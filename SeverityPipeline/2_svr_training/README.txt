# 2 — SVR Training

## Purpose

This part teaches a model to **predict the fault severity (in %)** from the features, then turns that continuous prediction into **severity classes**, evaluates everything honestly and saves a **model bundle** for later use.

## What goes in and what comes out

| **Input** | The feature database (optionally filtered), with feature names and severity labels |
| **Output** | A model bundle (`.mat`) + printed results + diagnostic figures |
| **Run** | `main_trainSvr.m` |

## Why Support Vector Regression?

SVR draws a smooth curve through the data. It keeps only the training examples that really matter (the *support vectors*) and ignores errors smaller than a tolerance. It works well with few samples and many features.

A **kernel** lets it follow curved relationships: a new signal is compared with the support vectors, and the closer it is to some of them, the more they influence the prediction.

## Workflow

```
 feature database
        │
        ▼
 1. Split train / test (stratified by severity)
        ▼
 2. Collapse feature families  (many windows → a summary)
        ▼
 3. Remove useless columns     (no variation)
        ▼
 4. Scale features             (comparable ranges)
        ▼
 5. Search model settings      (cross-validation, limited complexity)
        ▼
 6. Out-of-fold predictions    (honest predictions on training data)
        ▼
 7. Train final model
        ▼
 8. Evaluate on the untouched test set
        ▼
 9. Define severity classes
        ▼
10. Diagnostics, figures, save bundle
```

## Step by step

| Step | What happens | Why |

| Split | Part of the data is locked away as a test set | To measure performance on unseen data |
| Family collapse | Each feature observed over several windows is summarised by its average, and also its spread when the windows differ | Fewer, more robust inputs |
| Variance filter | Columns that never change are dropped | They carry no information |
| Scaling | Each feature is centred and rescaled | The kernel compares distances, so ranges must be comparable |
| Sample weighting | Healthy and severe samples count more | The most important cases must be learned well |
| Settings search | Model settings are tried automatically and scored by cross-validation | Avoid manual guessing |
| Complexity limit | Models that keep too many support vectors are penalised | Simpler models generalise better and run faster |
| Out-of-fold predictions | Every training sample is predicted by a model that never saw it | Needed to choose class boundaries without cheating |

**Important rule: every preprocessing step is learned on the training data only.** The test set never influences any decision; otherwise the evaluation would be too optimistic (data leakage).

## From a number to a class

The model outputs a continuous severity. To obtain classes:

1. A healthy / faulty limit is fixed in advance.
2. The fault range is divided into ordered classes by **dynamic programming**: the algorithm looks for the boundaries that make the predicted class match the true class as often as possible, under rules (minimum and maximum class size, preferred class width).
3. Several numbers of classes are compared; a small penalty favours simpler solutions.
4. The chosen boundaries are **frozen** and applied once to the test set.

## How performance is measured

| Level | Metrics |

| Regression (the %) | Errors (RMSE, MAE…), R² |
| Healthy vs faulty | Accuracy, precision, recall, specificity, F1, MCC |
| Severity classes | Confusion matrix, per-class precision/recall/F1, accuracy, balanced accuracy |
| Stability | Bootstrap resampling to see how sensitive the boundaries are |

## Diagnostics (sanity checks)

| Check | Question answered |

| Shuffled-label test | Does the model fail when labels are meaningless? (If not: leakage) |
| Nearest-neighbour test | Are test samples near-duplicates of training samples? |
| Support-vector analysis | How hard is the model "straining" to fit? |
| Noise robustness | Do predictions stay stable when features are slightly disturbed? |
| Learning curve | Would more data help? |
| Calibration | Are predicted severities consistent with true ones? |
| Feature importance | Which features matter most? |

## The saved bundle

| Content | Used for |

| Trained SVR model | Predict severity |
| Feature-family collapser | Reproduce step 2 on new data |
| Variance filter | Reproduce step 3 |
| Scaler | Reproduce step 4 |
| Names of the features used in training | Pick the right columns from a new signal |
| Class boundaries and class names | Convert severity to a class |
| Configuration | Traceability |

## Files

| File | Role |
|----------------------------------------------------------------------------------------|--------------------------------------|
| `main_trainSvr.m`                                                                      | Entry point                            |
| `SvrConfig.m`                                                                          | All settings in one place                  |
| `SeverityExperiment.m`                                                                 | Orchestrator of the whole workflow         |
| `SeverityDataset.m`, `StratifiedSplitter.m`                                            | Data loading and splitting |
| `SampleWeighter.m`                                                                     | Sample weights |
| `ConstrainedSvrSearch.m`, `TrialLog.m`, `CrossValidation.m`, `NestedCrossValidator.m`  | Settings search and validation |
| `OutOfFoldPredictor.m`                                                                 | Honest predictions on training data |
| `BoundaryEvaluator.m`, `BoundaryOptimizerDP.m`, `BoundaryCrossValidator.m`, `BoundarySelector.m`, `BootstrapStability.m` | Class boundary search |
| `RegressionMetrics.m`, `BinaryHealthMetrics.m`, `ConfusionMetrics.m` | Metrics |
| `ModelDiagnostics.m` | Sanity checks |
| `ConsoleReport.m`, `ResultPlotter.m`, `ConfusionPlotter.m` | Printed results and figures |
| `ModelBundleWriter.m` | Saves the bundle |

Shared pieces (scaler, family collapser, variance filter, SVR wrapper) are in `../common/`.

## Extending

| I want to… | Do this |

| Add a diagnostic | Add a method in `ModelDiagnostics` and call it from its run list |
| Add a figure | Add a method in `ResultPlotter` |
| Change training behaviour | Edit `SvrConfig` |
| Try another model | Replace the wrapper in `common/model` while keeping its two functions (train, estimate) |