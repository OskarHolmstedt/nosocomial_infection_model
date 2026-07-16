# Data and artifact policy

This policy distinguishes authored work from generated output so that the
repository remains understandable without risking the loss of research results.

## Version by default

- Source code, configuration, tests, and Quarto source documents.
- Small synthetic inputs needed by examples and tests.
- Selected figures and tables that are discussed in the README, a report, a
  poster, or a manuscript.
- Public poster material of reasonable size.
- Documentation needed to explain provenance and regeneration.

Selected outputs should have descriptive names and live under `results/figures/`
or `results/tables/`. Whenever practical, the generating script, parameter set,
and source data should be versioned with them.

The current poster is small enough for ordinary Git and is stored under
`docs/poster/` with a stable filename and a compressed PNG preview for web pages.

## Do not version by default

- R, RStudio, Posit, and operating-system session files.
- Quarto caches and generated support directories.
- Debug plots, generic interactive graphics, and temporary exports.
- Regenerable RDS caches and partial computations.
- Rendered reports unless they are intentionally published.
- Confidential or redistribution-restricted data.

These files belong in ignored directories such as `intermediate/`, `scratch/`,
or `tmp/`. Ignored does not mean disposable: do not delete an existing artifact
until its provenance and regeneration status have been confirmed.

## Large or archival outputs

Ordinary Git is appropriate for small, stable artifacts. Large binaries,
frequently changing outputs, and complete release bundles should be published as
GitHub Release assets or deposited in an archival repository such as Zenodo.
The repository should contain a link, version or DOI, and enough metadata to
connect the archive to the generating commit.

## Promoting an output to a versioned result

Before committing a generated artifact:

1. Confirm that it contains no confidential, embargoed, or restricted material.
2. Give it a descriptive, stable filename.
3. Move it to the appropriate `results/` or `docs/` location.
4. Record the script, parameters, input data, and commit used to produce it.
5. Prefer a compact web format for previews and PDF for print-ready figures.
