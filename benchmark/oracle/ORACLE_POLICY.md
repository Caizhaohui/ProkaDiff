# Oracle Environment Policy

See [docs/oracle_policy.md](../../docs/oracle_policy.md) for the complete scientific and parity policy.

This directory records pinned oracle toolchain versions (`versions.tsv`), benchmark datasets (`datasets.tsv`), and ground-truth oracle outputs (`expected/`).

Rule: Never update tool versions in `versions.tsv` without running and documenting a regression benchmark across the full oracle dataset suite.
