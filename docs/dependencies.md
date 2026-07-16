# Dependencies

The complete Conda environment is defined in [`environment.yaml`](../environment.yaml).
Create it from the repository root with Mamba or Conda:

```sh
mamba env create --file environment.yaml
conda activate nosocomial-infection-model
Rscript scripts/install_dependencies.R
```

Replace `mamba` with `conda` in the environment commands if Mamba is not
installed.

To update an existing environment after the specification changes:

```sh
mamba env update --file environment.yaml --prune
conda activate nosocomial-infection-model
Rscript scripts/install_dependencies.R
```

## Required for the quick start

- R 4.3 or newer;
- `here`, `Matrix`, `igraph`, `epicontacts`, `visNetwork`, and `expm` for the
  model and transmission representation;
- `ggplot2`, `dplyr`, `tidyr`, `scales`, and `htmlwidgets` for analysis and the
  interactive reconstruction.

Conda installs all core packages except `epicontacts`, which is not published on
the configured Conda channels. The idempotent `install_dependencies.R` step
installs its stable release from the package's
[official R-universe repository](https://reconhub.r-universe.dev/epicontacts)
and verifies the complete package inventory.

## Reports and diagnostics

- `patchwork` and `ggrepel` support diagnostic figures;
- `knitr`, `rmarkdown`, and Quarto render the reports;
- `testthat` runs the validation suite.

These are included in the complete environment because reports and validation
are part of the normal repository workflow, although the basic simulation and
inference functions do not require all of them.

## Optional high-resolution export

`webshot2`, `chromote`, `magick`, and `jsonlite` are used by
`save_timeline()` to export interactive timelines as print-ready PNG or PDF
files. They are included in the Conda environment, but export additionally needs
a Chrome or Chromium browser that `chromote` can locate.

The canonical quick start writes an HTML widget and does not require a browser
or the high-resolution export stack.
