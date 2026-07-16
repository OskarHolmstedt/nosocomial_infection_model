# Outbreak reconstruction for nosocomial infections

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Status: active research](https://img.shields.io/badge/status-active%20research-orange)

Research code for simulating hospital outbreaks and reconstructing plausible
transmission chains from incomplete epidemiological and genetic observations.

## Overview

Nosocomial outbreak investigations rarely observe every infected patient or the
exact time and source of each transmission. This project develops a Bayesian
framework for reasoning about those missing events. It combines a stochastic,
spatially structured hospital model with Markov chain Monte Carlo (MCMC)
inference of infection times, transmission ancestors, and unobserved generations
between detected cases.

The current work focuses on synthetic outbreaks, where the inferred transmission
history can be compared with known ground truth. The longer-term aim is to apply
the method to real hospital data.

## Model workflow

1. Construct a hospital with beds grouped into rooms and wards.
2. Simulate transmission, testing, admission, and discharge over time.
3. Retain the partially observed outbreak and pairwise genetic distances.
4. Use MCMC to sample compatible infection times and transmission trees.
5. Evaluate reconstruction quality against the simulated ground truth.

The model allows transmission rates to differ within a room, within a ward, and
across the hospital. Genetic information can be included or omitted to assess
how much it contributes to reconstruction.

## Example reconstruction

![Example posterior-mode reconstruction of a simulated hospital outbreak](code/figures/demo_mode.png)

Patient stays are arranged along the time axis and coloured by ward. Nodes mark
detected cases, while arrows show reconstructed transmission links. Edge colour
indicates the spatial scale of transmission. This example was generated from a
seeded synthetic outbreak.

## Repository guide

| Path | Contents |
| --- | --- |
| [`code/outbreak_simulation.R`](code/outbreak_simulation.R) | Stochastic outbreak simulation |
| [`code/mcmc.R`](code/mcmc.R) | Bayesian transmission reconstruction and scoring |
| [`code/prep_functions.R`](code/prep_functions.R) | Hospital layout and contact structure |
| [`code/parameters.R`](code/parameters.R) | Model parameters |
| [`code/main.R`](code/main.R) | Current simulation entry point |
| [`code/test.R`](code/test.R) | Current model checks |
| [`code/small_ward_data.R`](code/small_ward_data.R) | Example hospital configuration |
| [`environment.yaml`](environment.yaml) | Preliminary Conda environment specification |

Generated figures and experiment outputs are currently stored alongside the
research code. The repository layout and artifact policy are being reorganized in
[issue #4](https://github.com/OskarHolmstedt/nosocomial_infection_model/issues/4).

## Getting started

The current simulation entry point requires R with the `Matrix` and `igraph`
packages. Additional inference and reporting dependencies are being consolidated
as part of the reproducibility work.

Clone the repository, install the required packages, and run the current
demonstration from the `code` directory:

```sh
git clone https://github.com/OskarHolmstedt/nosocomial_infection_model.git
cd nosocomial_infection_model/code
Rscript main.R
```

This entry point simulates an outbreak and displays its infection-state matrix
and transmission tree. It is an exploratory workflow rather than a stable
command-line interface.

The committed Conda environment is not yet a complete dependency lockfile. A
smaller canonical quick start and a reproducible environment are tracked in
[issue #2](https://github.com/OskarHolmstedt/nosocomial_infection_model/issues/2)
and [issue #3](https://github.com/OskarHolmstedt/nosocomial_infection_model/issues/3).

## Project status

This repository contains research software under active development. Interfaces,
model assumptions, output formats, and results may change. The code has not yet
been prepared as a stable R package, and no archival software release has been
published.

No confidential hospital records or patient-level data should be committed to
this repository. The data currently included are generated from simulations.

## Project team

- Oskar Holmstedt — PhD student — [oskholms@chalmers.se](mailto:oskholms@chalmers.se)
- Philip Gerlee — supervisor — [gerlee@chalmers.se](mailto:gerlee@chalmers.se)
- Jon Edman Wallér — co-supervisor — [jon.edman@vgregion.se](mailto:jon.edman@vgregion.se)
- Torbjörn Lundh — co-supervisor — [torbjorn.lundh@chalmers.se](mailto:torbjorn.lundh@chalmers.se)

## Citation

There is not yet a citable software release or permanent DOI. Until one is
available, cite the repository and include the commit or release used:

> Holmstedt, O. *Outbreak reconstruction for nosocomial infections*. Research
> software, work in progress.
> <https://github.com/OskarHolmstedt/nosocomial_infection_model>

## License

The source code is available under the [MIT License](LICENSE).
