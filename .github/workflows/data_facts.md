# Data facts: OpenVaccine (stanford-covid-vaccine)

Source: the competition's **Data**, **Description** and **Evaluation** tabs (pasted text). Items marked **TO VERIFY** are not covered by that text and must be checked in the data itself.

## 1. Columns of `train.json` and what each means

| Column | Type | Meaning (from the Data page) | Role |
|---|---|---|---|
| `id` | string | Arbitrary identifier for each sample | key |
| `seq_scored` | integer | Number of positions with ground-truth values used in scoring. 68 in Train and Public Test, 91 in Private Test. Should match the length of `reactivity`, `deg_*` and `*_error_*` columns | metadata |
| `seq_length` | integer | Sequence length. 107 in Train and Public Test, 130 in Private Test | metadata |
| `sequence` | string | RNA sequence of A, G, U, C. 107 characters in Train and Public Test. First 68 bases correspond to the scored positions (indexed from 0) | input |
| `structure` | string | Dot-bracket notation of `(`, `)`, `.` describing whether each base is estimated to be paired or unpaired. Matching parentheses mark a base pair, e.g. `(....)` means base 0 pairs with base 5 and bases 1-4 are unpaired | input |
| `predicted_loop_type` | string | Structural context of each base, assigned by bpRNA from Vienna RNAfold 2 structure. S: paired Stem, M: Multiloop, I: Internal loop, B: Bulge, H: Hairpin loop, E: dangling End, X: eXternal loop | input |
| `reactivity` | list of floats, length 68 | Reactivity values for the first 68 bases; used to determine the likely secondary structure | **target (scored)** |
| `deg_pH10` | list of floats, length 68 | Degradation likelihood after incubating **without** magnesium at high pH (pH 10) | target (not scored) |
| `deg_Mg_pH10` | list of floats, length 68 | Degradation likelihood after incubating **with** magnesium at high pH (pH 10) | **target (scored)** |
| `deg_50C` | list of floats, length 68 | Degradation likelihood after incubating **without** magnesium at 50 degrees Celsius | target (not scored) |
| `deg_Mg_50C` | list of floats, length 68 | Degradation likelihood after incubating **with** magnesium at 50 degrees Celsius | **target (scored)** |
| `*_error_*` | list of floats, same length as the matching target | Calculated errors in the experimental values of `reactivity` and the `deg_*` columns | error estimates, not targets |
| S/N filter | flag | Indicates whether the sample passed the filters described in Section 5 | metadata |

**TO VERIFY in the data:** the exact name of the S/N filter column and any signal-to-noise value column (the page says "S/N filter" only). Run `df.columns.tolist()` and write down the real names. Also confirm the exact error column names (for example `reactivity_error`, `deg_error_Mg_pH10`).

## 2. Which target columns are scored?

The submission format requires all 5 targets, but **only three are scored**:

- `reactivity`
- `deg_Mg_pH10`
- `deg_Mg_50C`

Not scored: `deg_pH10` and `deg_50C` (the conditions without magnesium).

## 3. Which positions of each sequence are scored?

- **Train and Public Test:** the first **68** of 107 bases (positions 0-67). Measurements could not be carried out on the final bases for technical reasons.
- **Private Test:** the first **91** of 130 bases.
- Target arrays are therefore **shorter than the sequence** (68 vs 107), so labels must be aligned to the first 68 positions only.
- The sample submission lists a row for **every** position (`id_seqpos`, 457,953 rows), but only the scored positions count.

## 4. The metric (MCRMSE)

From the Evaluation tab: submissions are scored with **MCRMSE, mean columnwise root mean squared error**, where n_t is the number of scored ground-truth target columns, and y and y-hat are the actual and predicted values.

The formula image did not come through in the pasted text. In the standard form, consistent with the page's description:

    MCRMSE = (1 / n_t) * sum over scored columns j of  sqrt( (1 / n) * sum over scored positions i of (y_ij - yhat_ij)^2 )

- **n_t = 3**, the scored columns: `reactivity`, `deg_Mg_pH10`, `deg_Mg_50C`.
- Only positions up to `seq_scored` count. Positions beyond it are not scored.
- The two unscored columns (`deg_pH10`, `deg_50C`) must still be predicted for submissions but do not enter the metric.

**TO VERIFY:** compare the formula above with the rendered formula on the Evaluation tab, in particular that the inner sum runs over all scored positions of all sequences pooled (one RMSE per column), not per sequence.

Hand-worked example to turn into a unit test in Step 9 (fill in your own numbers):

> *(two sequences, three positions each, made-up predictions, RMSE per column, then the mean over the 3 columns)*

## 4b. Submission format

- One row per sequence position, with id `id_seqpos` (for example `id_00073f8be_0`), and a value for **all 5** targets: `reactivity, deg_Mg_pH10, deg_pH10, deg_Mg_50C, deg_50C`.
- A 107-long sequence needs 107 rows. Positions greater than `seq_scored` are not scored but still need a value.
- The sample submission has 457,953 rows (this equals 629 x 107 + 3005 x 130).
- AptaLens does not need to produce submission files, but keep the per-position layout in mind when you design the data model.

## 4c. Competition context

- Host: Stanford University (Rhiju Das, with the Eterna community). Goal: predict degradation rates at each base of an mRNA so that more stable vaccine candidates can be designed.
- Training data comes from a subset of an Eterna dataset of over 3000 RNA molecules. The final test sequences were synthesised and measured at Stanford in parallel to the competition.
- Timeline: entry and team-merger deadline 2 Oct 2020, final submission 6 Oct 2020 (extended by one day for a new data refresh). The competition is finished, so you can only make late submissions.
- Prizes: $25,000 total ($12,000, $8,000, $5,000 for the top three). 1,636 teams, 35,806 submissions.
- Tags: Biology, Biotechnology, Coronavirus, Public Health, Custom Metric.
- Citation: Das, Wayment-Steele, Kim, Choe, Tunguz, Reade, Demkin. *OpenVaccine: COVID-19 mRNA Vaccine Degradation Prediction.* Kaggle, 2020.

## 5. Noise filters and signal-to-noise

- **Data size:** 3029 sequences of length 107 originally. 629 were split out as the public test set and the remaining **2400** are in `train.json`.
- **Public test filtering.** The 629 public test sequences were filtered with three criteria:
  1. The minimum value across all 5 conditions must be greater than -0.5.
  2. Mean signal/noise across all 5 conditions must be greater than 1.0, where signal/noise = mean(measurement value over 68 nt) / mean(statistical error over 68 nt).
  3. For diversity, sequences were clustered at **less than 50% sequence similarity**, and test sequences were chosen from clusters with **3 or fewer members** (each test sequence is similar to at most 2 others).
- **These filters were NOT applied to `train.json`.** Training data can contain negative values, poor signal-to-noise, and **near-identical sequences**. The organisers provided everything in case competitors can extract more signal.
- **Private test:** about 3005 new RNAs of length 130 (measured on the first 91 bases). The page was updated on 18 Sep 2020: the three filters **will** be applied to the private test, and predictions on sequences that fail them are not scored.
- `test.json` contains 3634 sequences (629 public + 3005 private).

## 6. What this means for AptaLens

- **Train data has real leakage risk.** The organisers say it can contain near-identical sequences, which is the motivation for the similarity-aware splitters and leakage audit.
- **The competition's own public test set was built like a held-out-family split:** clusters at under 50% similarity, small clusters only. This supports using a ~50% identity threshold for the held-out-family splitter and an ~80% threshold for the clustered splitter.
- **Score filtering is a choice to make.** Decide whether to evaluate on all training sequences or only those passing the S/N and minimum-value filters, and report both if possible.
- **Distribution shift exists** between train (107 long, 68 scored) and private test (130 long, 91 scored). Phase 1 uses `train.json` only, so note this as a limitation.
- Always split by sequence `id`, never by row or position.

## 7. Facts to check against the data (Step 2 notebook)

- [ ] Number of rows is 2400 and `id` is unique
- [ ] `seq_length` is 107 and `seq_scored` is 68 for every row
- [ ] Target list lengths equal `seq_scored`
- [ ] Alphabets: sequence {A, C, G, U}, structure {`.`, `(`, `)`}, loop types {S, M, I, B, H, E, X}
- [ ] Target ranges and negative values (filter 1 was not applied to train)
- [ ] Names and distribution of the S/N filter and signal-to-noise columns
- [ ] Count of duplicate sequences, and near-duplicates (for example, shared first 30 bases)
- [ ] Formula on the Evaluation tab matches the one in section 4
- [ ] Hand-worked MCRMSE example written and checked
- [ ] Anything that differs from the description above