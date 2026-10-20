# AptaLens

A leakage-aware benchmarking and baseline toolkit for nucleic-acid sequence models, designed to match the scikit-learn-style API of pyaptamer.

## Why this exists

Reported scores for sequence models often collapse on sequences unlike the training data, because random train/validation splits leak near-identical sequences. AptaLens makes that visible and measurable:

- similarity-aware data splits,
- a leakage audit for any split,
- reproducible baselines and evaluation with confidence intervals.

## Who it's for

- Computational and wet-lab biologists screening aptamer or RNA candidates who need to know whether a model's score holds on novel sequences.
- pyaptamer contributors who want reusable splitters, metrics and benchmark runners.
- Anyone learning to build trustworthy ML evaluation for biological sequences.
