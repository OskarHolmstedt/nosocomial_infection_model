library(shiny)
library(bslib)
library(bsicons)
library(visNetwork)
library(here)

source(here("scripts", "_load_project.R"))

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "Nosocomial Outbreak Explorer",
  theme = bs_theme(preset = "shiny", version = 5),

  sidebar = sidebar(
    width = 300,

    accordion(
      open = c("Transmission", "Surveillance"),

      accordion_panel(
        "Transmission",
        icon = bs_icon("virus"),
        sliderInput("beta_room",
          HTML("&beta;<sub>room</sub> &mdash; same-room, per day"),
          min = 0.01, max = 0.30, value = 0.08, step = 0.01
        ),
        sliderInput("beta_ward",
          HTML("&beta;<sub>ward</sub> &mdash; same-ward, per day"),
          min = 0.000, max = 0.050, value = 0.005, step = 0.001
        ),
        sliderInput("beta_hosp",
          HTML("&beta;<sub>hosp</sub> &mdash; hospital-wide, per day"),
          min = 0.0000, max = 0.0050, value = 0.0005, step = 0.0001
        ),
        sliderInput("sim_days", "Simulation length (days)",
          min = 20, max = 90, value = 45, step = 5
        )
      ),

      accordion_panel(
        "Surveillance",
        icon = bs_icon("clipboard2-pulse"),
        sliderInput("test_rate", "Daily testing probability",
          min = 0.01, max = 0.50, value = 0.10, step = 0.01
        ),
        sliderInput("dis_rate", "Daily discharge probability",
          min = 0.05, max = 0.40, value = 0.10, step = 0.01
        )
      ),

      accordion_panel(
        "Inference",
        icon = bs_icon("cpu"),
        selectInput("n_samples", "MCMC samples",
          choices = c("500" = 500L, "1 000" = 1000L, "2 000" = 2000L),
          selected = 1000L
        ),
        numericInput("seed", "Random seed",
          value = 81369L, min = 1L, max = 99999L, step = 1L
        )
      )
    ),

    hr(),
    div(
      class = "d-grid gap-2",
      input_task_button("run", "Run simulation",
                        icon = bs_icon("play-fill"),
                        label_busy = "Running…",
                        class = "btn-primary")
    ),
    uiOutput("warn_msg")
  ),

  # ── Key metrics (live R0, post-run counts and accuracy) ────────────────────
  layout_column_wrap(
    width = 1 / 4, fill = FALSE,
    value_box(
      title   = "Basic reproduction number",
      value   = textOutput("r0_val", inline = TRUE),
      showcase = bs_icon("graph-up-arrow"),
      theme   = "primary",
      tooltip(
        bs_icon("info-circle", title = "About R\u2080"),
        "Updates live as you move the sliders. R\u2080 > 1 sustains an outbreak."
      )
    ),
    value_box(
      title   = "Total infections",
      value   = textOutput("total_cases", inline = TRUE),
      showcase = bs_icon("person-fill-exclamation"),
      theme   = "secondary"
    ),
    value_box(
      title   = "Observed (test +ve)",
      value   = textOutput("obs_cases", inline = TRUE),
      showcase = bs_icon("eyedropper"),
      theme   = "info"
    ),
    value_box(
      title   = "Mode accuracy",
      value   = textOutput("accuracy", inline = TRUE),
      showcase = bs_icon("bullseye"),
      theme   = "success",
      tooltip(
        bs_icon("info-circle", title = "About mode accuracy"),
        "Fraction of cases where the posterior-mode ancestor matches the true ancestor."
      )
    )
  ),

  # ── Timeline cards ─────────────────────────────────────────────────────────
  layout_column_wrap(
    width = 1 / 2,

    card(
      full_screen = TRUE,
      card_header("Ground-truth transmission tree"),
      visNetworkOutput("true_tree", height = "430px")
    ),

    card(
      full_screen = TRUE,
      card_header(
        "Posterior-mode reconstruction",
        tooltip(
          bs_icon("info-circle", title = "Reading the reconstruction"),
          HTML("Green edges = correctly inferred.<br>
                Red edges = errors.<br>
                Edge width \u221d posterior support.")
        )
      ),
      visNetworkOutput("recon_tree", height = "430px")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive copies of simulation results
  outbreak_rv <- reactiveVal(NULL)
  fit_rv      <- reactiveVal(NULL)
  score_rv    <- reactiveVal(NULL)
  warn_rv     <- reactiveVal(NULL)

  # Build a modified theta from slider values (no simulation involved)
  build_theta <- function() {
    th            <- theta                          # global from _load_project.R
    th$beta["BrB"] <- input$beta_room
    th$beta["BwB"] <- input$beta_ward
    th$beta["BhB"] <- input$beta_hosp
    th$P$Test      <- rep(input$test_rate, NumBeds)
    th$P$Dis       <- rep(input$dis_rate,  NumBeds)
    th$P$Init      <- c(1, rep(0, NumBeds - 1), rep(0, NumRooms))
    W              <- weights(th$beta, Contact)
    th$LogComp     <- log1p(-W)
    th$Odds        <- W / (1 - W)
    th$spatial_sizes <- theta$spatial_sizes         # unchanged
    th
  }

  # R0 updates live without running the simulation
  output$r0_val <- renderText({
    req(input$beta_room, input$beta_ward, input$beta_hosp, input$dis_rate)
    th <- build_theta()
    sprintf("%.2f", r_0(th))
  })

  # Placeholder text before first run
  output$total_cases <- renderText({ if (is.null(outbreak_rv())) "—" else nrow(outbreak_rv()$CaseRec) })
  output$obs_cases   <- renderText({ if (is.null(outbreak_rv())) "—" else nrow(outbreak_rv()$ObsRec) })
  output$accuracy    <- renderText({
    sc <- score_rv()
    if (is.null(sc)) "—" else sprintf("%.0f%%", sc$mode_accuracy * 100)
  })

  # ── Main pipeline triggered by the Run button ──────────────────────────────
  observeEvent(input$run, {
    warn_rv(NULL)

    th <- build_theta()

    # 1 — Simulate
    seed <- as.integer(isolate(input$seed))
    set.seed(seed, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    ob <- simulate_outbreak(as.integer(isolate(input$sim_days)), th)
    outbreak_rv(ob)

    # Guard: need at least 2 observed cases to run MCMC
    if (nrow(ob$ObsRec) < 2L) {
      warn_rv("Outbreak died out — fewer than 2 cases detected. Try higher \u03b2 or lower discharge rate.")
      fit_rv(NULL)
      score_rv(NULL)
      return()
    }

    # 2 — MCMC
    n  <- as.integer(isolate(input$n_samples))
    bi <- max(100L, as.integer(n * 0.2))
    set.seed(seed + 1L, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    fit <- mcmc(ob, N_samples = n, burn_in = bi)
    fit_rv(fit)

    # 3 — Score
    sc <- ancestry_score(ob$ObsRec$Anc2, fit$anc,
                         ob$ObsRec$Adm,  ob$ObsRec$PTest)
    score_rv(sc)

    show_toast(
      toast(
        sprintf("Done — %d cases, %.0f%% accuracy",
                nrow(ob$ObsRec), sc$mode_accuracy * 100),
        header = "Simulation complete",
        type   = "success"
      )
    )
  })

  # ── Warning message ────────────────────────────────────────────────────────
  output$warn_msg <- renderUI({
    msg <- warn_rv()
    if (is.null(msg)) return(NULL)
    div(class = "mt-2 text-warning small", bs_icon("exclamation-triangle"), msg)
  })

  # ── Timeline outputs ───────────────────────────────────────────────────────
  output$true_tree <- renderVisNetwork({
    ob <- outbreak_rv()
    req(ob, nrow(ob$ObsRec) >= 2L)
    plot_timeline(ob, observed = TRUE, show_all_stays = TRUE)
  })

  output$recon_tree <- renderVisNetwork({
    ob  <- outbreak_rv()
    fit <- fit_rv()
    req(ob, fit)
    plot_timeline_mode(ob, fit, color_edges = TRUE, color_nodes = FALSE)
  })
}

shinyApp(ui, server)
