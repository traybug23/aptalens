# The biology behind AptaLens

AptaLens is a toolkit for evaluating machine-learning models that predict properties of nucleic-acid sequences. This document explains the biology the project rests on: what aptamers are, which properties are worth predicting, what the models are scored on, and why the evaluation is designed the way it is.

## 1. Motivation

Aptamers are **short single-stranded DNA or RNA molecules** that fold into three-dimensional shapes and bind specific targets, such as proteins or small molecules. They are found by iterative laboratory selection and, increasingly, designed or screened computationally.

Computational aptamer work depends on models that map **sequence to property**. 

Those models are only useful if their reported performance holds on sequences unlike the ones they were trained on. 

In aptamer data, that is easy to get wrong: selection experiments produce large pools of closely related sequences, so a random train/validation split puts near-duplicates on both sides and inflates scores.

AptaLens exists to make that problem measurable. It provides similarity-aware data splits, a leakage audit, reproducible baselines, and evaluation with confidence intervals, with an API like pyaptamer.

## 2. Aptamer biology in brief

- **Structure determines function.** An aptamer's sequence folds into stems (base-paired helices) and loops (unpaired regions), which assemble into a 3D shape. Binding depends on that shape, so sequence-to-structure and structure-to-function are the core relationships any model has to learn.
- **How aptamers are found.** In SELEX (Systematic Evolution of Ligands by EXponential enrichment), a large random pool of sequences is repeatedly exposed to a target, binders are kept and amplified, and the cycle is repeated. The result is a set of related sequences enriched over rounds, which is why **sequence families and near-duplicates are the norm, not the exception**.
- **Stability is a practical requirement.** Nucleic acids are degraded by chemical hydrolysis and, in biological fluids, by nucleases. RNA is especially vulnerable because each sugar carries a 2'-hydroxyl group that can attack the neighbouring backbone phosphate and cut the chain. Chemical modifications are widely used to improve stability. An aptamer that cannot survive its environment is of little use, however well it binds.

## 3. Why position-level stability and accessibility

Binding affinity is the property most people want, but good affinity data is scarce and noisy. Position-level measurements of **flexibility** and **degradation** are different in three ways:

1. They describe the molecule's physical behaviour at every nucleotide, which links directly to structure: bases locked in stems are rigid and protected, while unpaired bases in loops and ends are more flexible and exposed.
2. They are available in larger, more consistent datasets, which makes them a better place to develop and test evaluation methodology.
3. They are practically relevant, because stability is a real design constraint, and the same sequence-to-property modeling skills transfer to binding prediction.

## 4. What AptaLens predicts, and what is scored

**Input per molecule:** a nucleotide sequence, plus its predicted secondary structure and per-base loop type (stem, hairpin, internal loop, bulge, multiloop, external loop, dangling end).

**Output per nucleotide:** three measured quantities.

| Target | What it reflects biologically |
| --- | --- |
| `reactivity` | How flexible or accessible a nucleotide is. High values point to unpaired, exposed positions, so this is a readout of structure. |
| `deg_Mg_pH10` | Degradation tendency at each position at high pH in the presence of magnesium. |
| `deg_Mg_50C` | Degradation tendency at each position at elevated temperature in the presence of magnesium. |

Two further measurements (the same degradation conditions without magnesium) are recorded in the dataset but are not part of the score.

**Scoring.** Models are scored with **MCRMSE**, the mean over the three targets of the root mean squared error per target, computed over the positions that have measurements. Each target also has a measurement-error estimate, which is why uncertainty is reported alongside scores.

![MCRMSE formula]("C:\Users\sneha\OneDrive\Pictures\Screenshots\Screenshot 2026-10-23 005717.png")

**Evaluation protocol.** Each model is scored under three split strategies: random, similarity-clustered, and held-out family. The drop in score from the first to the last is the measure of how much a model relies on memorizing close relatives instead of learning generalizable sequence-structure rules.

## 5. Why the evaluation is designed this way

- **Splits are made by molecule, never by position.** Positions within one molecule are strongly dependent, so splitting them apart leaks information.
- **Similarity-aware clustering.** Sequences are grouped by sequence identity, and whole clusters are assigned to either train or validation. A stricter identity threshold (around 50%) approximates "unrelated families", and a looser one (around 80%) removes near-duplicates.
- **Leakage audit.** For each validation sequence, AptaLens computes the identity to its closest training sequence and reports the distribution. This turns "is this split fair?" from an assumption into a number.
- **Bootstrap confidence intervals** resample whole molecules, not positions, because positions are not independent samples.
- **Error analysis uses biology.** Breakdowns by loop type, position along the molecule, GC content, and measurement noise show where and why a model fails, for example whether flexible loop regions are harder than stems.

## 6. Data used in Phase 1

Phase 1 uses a public dataset of 2,400 RNA molecules of length 107, with measurements at the first 68 positions under five conditions, of which three are scored. It was released with the Stanford / Eterna OpenVaccine challenge (Kaggle, `stanford-covid-vaccine`). It was chosen as a proxy task because it is large, clean, per-nucleotide, and rooted in the same structure-stability relationships described above. The project's design does not depend on it: a dataset registry lets aptamer datasets be added, and moving the pipeline to aptamer data is the goal of Phase 3.

Properties of this data that matter for the biology:

- The structure and loop-type features are **computationally predicted**, not experimentally determined, so models inherit their errors.
- The sequences are human-designed, so many are close variants of each other. This is the same near-duplicate structure found in selection pools, and a good test of the leakage tools.
- Measurements exist only for the first 68 of 107 positions, so labels must be aligned to those positions only.

## 7. Scope and limitations

- This is a **proxy task** (stability and accessibility), not binding affinity. Performance here does not imply performance on binding.
- The proxy data is **RNA**, while many therapeutic aptamers are DNA or chemically modified. Transfer between them is an open question this project does not claim to settle.
- Per-position degradation depends on experimental conditions (pH, temperature, ions) that differ from physiological ones.
- Predicted structures are imperfect, and real aptamers can adopt more than one fold.
- Small, noisy aptamer datasets limit what can be concluded, which is why the project reports confidence intervals and states its assumptions.

## 8. Glossary

- **Aptamer:** a short single-stranded nucleic acid that folds to bind a specific target.
- **SELEX:** laboratory selection process used to find aptamers.
- **Secondary structure:** the pattern of base pairing within a strand.
- **Stem / loop:** a paired helical region / an unpaired region.
- **Dot-bracket notation:** text encoding of pairing, where `.` is unpaired and matching `(` and `)` are paired bases.
- **Reactivity:** a measure of a base's flexibility or accessibility.
- **Degradation:** breakage of the nucleic-acid backbone.
- **MCRMSE:** mean columnwise root mean squared error.
- **Leakage:** information from the validation set reaching training, usually through near-duplicate sequences.
- **Held-out family:** a split in which entire clusters of related sequences are reserved for validation.