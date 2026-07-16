# Data

This directory is reserved for documented inputs that are safe and useful to
share with the repository.

## What may be committed

- Small simulated datasets used by examples or tests.
- Stable, processed research outputs when they are needed to reproduce a
  versioned result and are not practical to regenerate during a normal run.
- Metadata and data dictionaries that do not disclose confidential information.

Every committed dataset should state how it was created, which script and
parameters produced it, and the date or commit associated with it. Prefer open,
portable formats such as CSV for tabular data. An RDS file may accompany a CSV
when it materially improves the R workflow, but should not be the only copy of
important tabular results.

## What must not be committed

- Patient-level hospital records or other confidential data.
- Direct or indirect identifiers.
- Data governed by an agreement that does not permit public redistribution.
- Credentials, access tokens, or local paths to restricted storage.

Restricted data should remain in an approved external location. Code should
accept its location through configuration outside version control; do not encode
the storage path or access details in the repository.

## Generated and intermediate data

Regenerable caches and partial computations belong in `intermediate/`, which is
ignored by Git. Selected aggregate results intended for presentation belong in
`results/` and follow the repository's
[artifact policy](../docs/artifact-policy.md).
