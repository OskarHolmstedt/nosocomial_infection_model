library(shiny)
library(bslib)
library(bsicons)
library(visNetwork)
library(ggplot2)
library(here)

source(here("scripts", "_load_project.R"))

# ── Helpers ───────────────────────────────────────────────────────────────────

num_slider <- function(id, label, min, max, value, step, ...) {
  sliderInput(id, label, min = min, max = max, value = value,
              step = step, ticks = FALSE, ...)
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "Nosocomial Outbreak Explorer",
  theme = bs_theme(preset = "shiny", version = 5),

  sidebar = sidebar(
    width = 320,
    style = "overflow-y: auto; max-height: calc(100vh - 56px);",

    accordion(
      open = c("Patient transmission", "Surveillance"),

      # ── Hospital layout ──────────────────────────────────────────────────
      accordion_panel(
        "Hospital layout",
        icon = bs_icon("building"),
        num_slider("wards",          "Wards",          1, 4, 2, 1),
        num_slider("rooms_per_ward", "Rooms per ward",  1, 8, 5, 1),
        num_slider("beds_per_room",  "Beds per room",   2, 6, 4, 1),
        helpText("Changes here rebuild the contact matrix before each run.")
      ),

      # ── Patient-to-patient transmission ──────────────────────────────────
      accordion_panel(
        "Patient transmission",
        icon = bs_icon("arrow-left-right"),
        num_slider("beta_room", HTML("&beta;<sub>room</sub> — same room, per day"),
                   0.01, 0.30, 0.08, 0.01),
        num_slider("beta_ward", HTML("&beta;<sub>ward</sub> — same ward, per day"),
                   0.000, 0.050, 0.005, 0.001),
        num_slider("beta_hosp", HTML("&beta;<sub>hosp</sub> — hospital-wide, per day"),
                   0.0000, 0.0050, 0.0005, 0.0001)
      ),

      # ── Room transmission (feature-gated) ────────────────────────────────
      accordion_panel(
        "Room transmission",
        icon = bs_icon("door-open"),
        input_switch("use_room_transmission",
                     "Enable patient ↔ room spread", value = FALSE),
        conditionalPanel(
          "input.use_room_transmission",
          num_slider("beta_br", HTML("&beta;<sub>BR</sub> — patient → room, per day"),
                     0.00, 0.20, 0.00, 0.01),
          num_slider("beta_rb", HTML("&beta;<sub>RB</sub> — room → patient, per day"),
                     0.00, 0.20, 0.00, 0.01),
          input_switch("use_room_room",
                       "Enable room → room spread", value = FALSE),
          conditionalPanel(
            "input.use_room_transmission && input.use_room_room",
            num_slider("beta_rr", HTML("&beta;<sub>RR</sub> — room → room, per day"),
                       0.00, 0.10, 0.00, 0.01)
          )
        )
      ),

      # ── Surveillance ──────────────────────────────────────────────────────
      accordion_panel(
        "Surveillance",
        icon = bs_icon("clipboard2-pulse"),
        num_slider("test_rate",     "Daily patient testing prob.",     0.01, 0.50, 0.10, 0.01),
        num_slider("dis_rate",      "Daily discharge prob.",           0.05, 0.40, 0.10, 0.01),
        num_slider("clean_rate",    "Daily decontamination prob.",     0.01, 0.30, 0.05, 0.01),
        num_slider("room_test_rate","Daily room swab prob.",           0.00, 0.20, 0.05, 0.01),
        num_slider("imp_rate",      "Daily community importation prob.",0.00, 0.20, 0.00, 0.01)
      ),

      # ── Initial conditions ────────────────────────────────────────────────
      accordion_panel(
        "Initial conditions",
        icon = bs_icon("play-circle"),
        num_slider("n_seeds",  "Seed infections at t = 0", 1, 10, 1, 1),
        num_slider("sim_days", "Simulation length (days)", 20, 90, 45, 5)
      ),

      # ── Genomics ──────────────────────────────────────────────────────────
      accordion_panel(
        "Genomics",
        icon = bs_icon("braces-asterisk"),
        num_slider("mu", HTML("&mu; — mutations per transmission"), 1, 10, 2, 1)
      ),

      # ── MCMC priors ───────────────────────────────────────────────────────
      accordion_panel(
        "MCMC priors",
        icon = bs_icon("diagram-3"),
        num_slider("lambda_import",
                   HTML("&lambda;<sub>import</sub> — geometric import rate"),
                   0.01, 0.60, 0.30, 0.01),
        num_slider("pi_import",
                   HTML("&pi;<sub>import</sub> — community import probability"),
                   0.00, 0.30, 0.05, 0.01),
        helpText(HTML("&lambda;<sub>import</sub> should be &ge; testing rate
                       to overcome the testing-likelihood pull on import times."))
      ),

      # ── Inference ─────────────────────────────────────────────────────────
      accordion_panel(
        "Inference",
        icon = bs_icon("cpu"),
        selectInput("n_samples", "MCMC samples",
                    choices = c("500" = 500L, "1 000" = 1000L,
                                "2 000" = 2000L, "3 000" = 3000L),
                    selected = 1000L),
        input_switch("use_rooms_mcmc",
                     "Room-aware inference (patients_and_rooms)", value = FALSE),
        numericInput("seed", "Random seed",
                     value = 81369L, min = 1L, max = 99999L, step = 1L)
      )
    ),

    hr(),
    div(
      class = "d-grid gap-2",
      input_task_button("run", "Run simulation",
                        icon         = bs_icon("play-fill"),
                        label_busy   = "Running…",
                        class        = "btn-primary")
    ),
    uiOutput("warn_msg")
  ),

  # ── Key metrics (R0 updates live; others after Run) ───────────────────────
  layout_column_wrap(
    width = 1 / 4, fill = FALSE,
    value_box(
      title    = "Basic reproduction number",
      value    = textOutput("r0_val", inline = TRUE),
      showcase = bs_icon("graph-up-arrow"),
      theme    = "primary",
      tooltip(
        bs_icon("info-circle", title = "About R\u2080"),
        "Updates live. Values > 1 sustain an outbreak."
      )
    ),
    value_box(
      title    = "Total infections",
      value    = textOutput("total_cases", inline = TRUE),
      showcase = bs_icon("person-fill-exclamation"),
      theme    = "secondary"
    ),
    value_box(
      title    = "Observed (test +ve)",
      value    = textOutput("obs_cases", inline = TRUE),
      showcase = bs_icon("eyedropper"),
      theme    = "info"
    ),
    value_box(
      title    = "Mode accuracy",
      value    = textOutput("accuracy", inline = TRUE),
      showcase = bs_icon("bullseye"),
      theme    = "success",
      tooltip(
        bs_icon("info-circle", title = "About mode accuracy"),
        "Fraction of cases where the posterior-mode ancestor matches the true ancestor."
      )
    )
  ),

  # ── Output cards ──────────────────────────────────────────────────────────
  layout_column_wrap(
    width = 1 / 3,

    card(
      full_screen = TRUE,
      card_header("Prevalence dynamics"),
      plotOutput("prev_plot", height = "350px")
    ),

    card(
      full_screen = TRUE,
      card_header("Ground-truth transmission tree"),
      visNetworkOutput("true_tree", height = "350px")
    ),

    card(
      full_screen = TRUE,
      card_header(
        "Posterior-mode reconstruction",
        tooltip(
          bs_icon("info-circle", title = "Reading the reconstruction"),
          HTML("Green = correctly inferred.<br>Red = error.<br>
                Width \u221d posterior support.")
        )
      ),
      visNetworkOutput("recon_tree", height = "350px")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  outbreak_rv <- reactiveVal(NULL)
  fit_rv      <- reactiveVal(NULL)
  score_rv    <- reactiveVal(NULL)
  warn_rv     <- reactiveVal(NULL)

  # ── Hospital layout (reactive, rebuilt when layout sliders change) ─────────
  hosp_rv <- reactive({
    hospital(
      Wards        = as.integer(input$wards),
      RoomsPerWard = as.integer(input$rooms_per_ward),
      BedsPerRoom  = as.integer(input$beds_per_room)
    )
  })

  # ── Build theta from all slider values ────────────────────────────────────
  build_theta <- function(h) {
    nb <- h$NumBeds
    nr <- h$NumRooms

    use_room <- isTRUE(input$use_room_transmission)

    beta <- c(
      BrB = input$beta_room,
      BwB = input$beta_ward,
      BhB = input$beta_hosp,
      BR  = if (use_room) input$beta_br else 0,
      RB  = if (use_room) input$beta_rb else 0,
      RR  = if (use_room && isTRUE(input$use_room_room)) input$beta_rr else 0
    )

    W <- weights(beta, h$Contact)

    n_seeds <- min(as.integer(input$n_seeds), nb)
    init    <- c(rep(1, n_seeds), rep(0, nb - n_seeds), rep(0, nr))

    list(
      mu     = as.integer(input$mu),
      beta   = beta,
      P = list(
        Init     = init,
        Imp      = rep(input$imp_rate,       nb),
        Test     = rep(input$test_rate,      nb),
        Dis      = rep(input$dis_rate,       nb),
        Clean    = rep(input$clean_rate,     nr),
        RoomTest = rep(input$room_test_rate, nr)
      ),
      use_room_transmission = use_room,
      use_room_room         = use_room && isTRUE(input$use_room_room),
      lambda_import         = input$lambda_import,
      pi_import             = input$pi_import,
      LogComp               = log1p(-W),
      Odds                  = W / (1 - W),
      spatial_sizes         = h$spatial_sizes,
      Positions             = h$Positions
    )
  }

  # ── R0: live, no run needed ───────────────────────────────────────────────
  output$r0_val <- renderText({
    req(input$beta_room, input$beta_ward, input$beta_hosp, input$dis_rate)
    th <- build_theta(hosp_rv())
    sprintf("%.2f", r_0(th))
  })

  # Placeholders before first run
  output$total_cases <- renderText({
    ob <- outbreak_rv(); if (is.null(ob)) "—" else nrow(ob$CaseRec)
  })
  output$obs_cases <- renderText({
    ob <- outbreak_rv(); if (is.null(ob)) "—" else nrow(ob$ObsRec)
  })
  output$accuracy <- renderText({
    sc <- score_rv()
    if (is.null(sc)) "—" else sprintf("%.0f%%", sc$mode_accuracy * 100)
  })

  # ── Main pipeline ─────────────────────────────────────────────────────────
  observeEvent(input$run, {
    warn_rv(NULL)

    h  <- isolate(hosp_rv())
    th <- isolate(build_theta(h))

    # 1 ── Simulate
    seed <- as.integer(isolate(input$seed))
    set.seed(seed, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    ob <- simulate_outbreak(as.integer(isolate(input$sim_days)), th)
    outbreak_rv(ob)

    if (nrow(ob$ObsRec) < 2L) {
      warn_rv("Outbreak died out — fewer than 2 cases detected. Try raising \u03b2 or lowering discharge rate.")
      fit_rv(NULL); score_rv(NULL)
      return()
    }

    # 2 ── MCMC
    n    <- as.integer(isolate(input$n_samples))
    bi   <- max(100L, as.integer(n * 0.2))
    mode <- if (isolate(input$use_rooms_mcmc) && isTRUE(th$use_room_transmission))
              "patients_and_rooms" else "patients"

    set.seed(seed + 1L, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    # use_room_tests = FALSE: all simulated room episodes are candidates.
    # TRUE would restrict to swab-positive rooms only; at 5% daily swab rate
    # most room ancestors would be invisible and accuracy collapses. Both the
    # MCMC call and the scoring call must use the same setting.
    use_rt <- FALSE
    fit <- mcmc(ob, N_samples = n, burn_in = bi,
                inference_mode = mode, use_room_tests = use_rt)
    fit_rv(fit)

    # 3 ── Score (patient-only or room-aware depending on inference mode)
    sc <- if (mode == "patients_and_rooms") {
      truth <- build_room_aware_truth(ob, use_room_tests = use_rt)
      ancestry_score(truth$true_anc, fit$anc,
                     truth$adm_times, truth$ptest_times)
    } else {
      ancestry_score(ob$ObsRec$Anc2, fit$anc,
                     ob$ObsRec$Adm,  ob$ObsRec$PTest)
    }
    score_rv(sc)

    show_toast(
      toast(
        sprintf("%d cases detected, %.0f%% reconstruction accuracy",
                nrow(ob$ObsRec), sc$mode_accuracy * 100),
        header = "Simulation complete",
        type   = "success"
      )
    )
  })

  # ── Warning banner ────────────────────────────────────────────────────────
  output$warn_msg <- renderUI({
    msg <- warn_rv()
    if (is.null(msg)) return(NULL)
    div(class = "mt-2 text-warning small",
        bs_icon("exclamation-triangle-fill"), " ", msg)
  })

  # ── Prevalence plot ───────────────────────────────────────────────────────
  output$prev_plot <- renderPlot({
    ob <- outbreak_rv(); req(ob)
    plot_prevalence(outbreak_prevalence(ob))
  }, res = 96)

  # ── Ground-truth timeline ─────────────────────────────────────────────────
  output$true_tree <- renderVisNetwork({
    ob <- outbreak_rv(); req(ob, nrow(ob$ObsRec) >= 2L)
    plot_timeline(ob, observed = TRUE, show_all_stays = TRUE,
                  show_rooms = isTRUE(input$use_room_transmission))
  })

  # ── Reconstruction timeline ───────────────────────────────────────────────
  output$recon_tree <- renderVisNetwork({
    ob  <- outbreak_rv(); req(ob)
    fit <- fit_rv();      req(fit)
    plot_timeline_mode(ob, fit, color_edges = TRUE, color_nodes = FALSE)
  })
}

shinyApp(ui, server)
