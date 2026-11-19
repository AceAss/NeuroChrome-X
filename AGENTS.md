# NeuroChrom-X: rules for coding agents
Read docs/PRD.md and the Task Card before acting. The PRD and the ADRs are the spec.

## Process
- Work only from a Task Card. One card = one branch. Small commits with conventional messages.
- Never change thresholds, expected values in tests, ADRs or the PRD. If a test fails, report it; do not edit the test to make it pass.
- Never edit files that a Task Card lists as provided by Claude. Report problems instead.
- Dry-run first (snakemake -n), then toy data, then real data.
- Use only documented tool flags. If unsure, run the tool with --help and quote the relevant lines in the report. Do not guess flags.
- Pin every dependency and update the environment lockfiles. List any new dependency in the report.

## Data and safety
- Never overwrite raw data. All outputs go under results/. No hard-coded paths; use config.
- Ask before downloading more than 20 GB or more than 50 files. Record checksums and provenance for every download.
- Stop and ask if the genome build or annotation is ambiguous, a QC gate fails, or a step needs more than about 4 hours or 32 GB of RAM (adjust after the machine audit).
- Never commit secrets. Read API keys (for example for AlphaGenome) from environment variables and keep them out of logs and reports.
- For external models or APIs, record the model or API version and the call date with every output, cache results, and stop if you hit rate limits. Do not start GPU or cloud jobs that cost money without approval.

## Code
- Python 3.11 or later (ruff, black, pytest) with type hints on public functions. R scripts for Bioconductor steps, called from Snakemake.

## Reporting
- Finish every task with a Run Report in the format given in the Task Card.
