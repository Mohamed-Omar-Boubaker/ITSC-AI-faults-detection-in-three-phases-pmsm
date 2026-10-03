# 1 — Database Generation

## Purpose

A machine-learning model cannot read raw currents directly. This part converts every recorded sample into a **row of numbers (features)** and stores all rows in a table: the **feature database**.

## What goes in and what comes out
|---------------------------------------------------------------------------------------------------|
|            |                                  Description                                         |
|------------|--------------------------------------------------------------------------------------|
| **Input**  | A file of raw three-phase current samples, each with its known fault severity        |
| **Output** | A `.mat` table: one row per sample, one column per feature, plus the severity labels |
| **Run**    | `main_buildDatabase.m`                                                               |

## How it works

```
 raw sample (ia, ib, ic) + severity label
          │
          ▼
   ┌─────────────┐
   │ Resampling  │  same rate for every sample
   └──────┬──────┘
          ▼
   ┌───────────────────────────────┐
   │ Two extraction pipelines      │
   │                               │
   │  A) Cyclic window features    │──┐
   │  B) Wavelet energy features   │──┤
   └───────────────────────────────┘  │
                                      ▼
                              merge into one row
                                      │
                                      ▼
                         saved after every sample
                         (a long run can be resumed)
```

## Key concepts

**Resampling.** Raw data may be recorded faster than needed. It is converted once to a common rate, so both pipelines work on the same signal.

**Symmetrical components.** The three currents are decomposed into a balanced part (positive sequence) and an imbalance part (negative sequence). The ratio between them is a strong fault indicator.

**Cyclic analysis.** Instead of one value per sample, every feature is computed on each consecutive window. The feature is then described by a short list of values: its evolution in time.

## The two pipelines

| Pipeline | Works on | Features produced |
|---|---|---|
| **A — Cyclic window features** | Positive-sequence, negative-sequence and space-vector signals, window by window | Time-domain, frequency-domain, geometric, nonlinear features and the current unbalance factor |
| **B — Wavelet energy features** | The space vector, using a bilateral wavelet transform with synchrosqueezing | Energy reconstructed around the fault frequency, wavelet entropy, sub-band energy ratio, all per window |

## Feature families

| Family | What it tells |
|---|---|
| Time-domain | Amplitude, spread, peakiness, smoothness of the waveform |
| Frequency-domain | Distribution of energy over frequencies, harmonic content |
| Geometric | Shape of the space-vector orbit (circle when balanced, ellipse when not) |
| Nonlinear | Irregularity and complexity of the waveform |
| Unbalance factor | How strong the imbalance is compared with the balanced part |
| Wavelet energy | How much energy sits in the band where the fault leaves its trace |

## Structure of the output table

```
             feature 1         feature 2        ...
           ┌────┬────┬────┐ ┌────┬────┬────┐
 sample 1  │ w1 │ w2 │ w3 │ │ w1 │ w2 │ w3 │ ...
 sample 2  │ ...                              
```

Each column name is *feature name* + *window index*. Columns sharing the same feature name form a **family**; the training part later summarises each family.

## Files
|--------------------------------------------------------------------------------------------------------|
|                File                        |                           Role                            |
|--------------------------------------------|-----------------------------------------------------------|
| `main_buildDatabase.m`                     | Entry point: assembles the parts and runs the build       |
| `FeatureDatabaseBuilder.m`                 | Orchestrator: read → resample → extract → merge → save    |
| `SampleSource.m` / `Hdf5SampleSource.m`    | Interface and implementation for reading raw samples      |
| `DatabaseWriter.m` / `MatDatabaseWriter.m` | Interface and implementation for saving the table         |

The feature extractors themselves live in `../common/features/`.

