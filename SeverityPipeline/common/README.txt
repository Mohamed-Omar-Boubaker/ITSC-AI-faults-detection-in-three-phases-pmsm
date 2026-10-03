# Common — Shared Building Blocks

## Purpose

Several parts of the project need exactly the same tools. Keeping them in one place guarantees that **the database, the training and the inference all treat signals identically**, and avoids duplicated code that could drift apart.

## Structure

```
common/
├── signal/          handling of currents and signals
├── features/        feature extraction
├── preprocessing/   preparation of feature tables
└── model/           the learning model wrapper
```

## Contents

| Subfolder | Contents | Used by |

| `signal/` | Phase currents container, resampling, symmetrical components, signal bundle, windowing | 1, 3 |
| `features/` | Feature groups (time, frequency, geometric, nonlinear), window-level blocks, the two extraction pipelines, wavelet and synchrosqueezing tools, band tracking, factory that defines the feature set | 1, 3 |
| `preprocessing/` | Feature-family collapser, variance filter, standard scaler | 2, 3 |
| `model/` | SVR wrapper (train and estimate) | 2, 3 |

## How the feature extraction is organised

```
FeatureExtractorFactory   ← the single place that decides WHICH features exist
        │
        ├── CyclicFeatureExtractor  (pipeline A)
        │       └── window blocks → feature groups (time, frequency, geometric, nonlinear, unbalance)
        │
        └── SwtFeatureExtractor     (pipeline B)
                └── wavelet transform → synchrosqueezing → band tracking → energy metrics
```

| Class type | Role |

| Feature group | A list of named features computed on one signal; adding a feature = adding a line |
| Window feature | Applies a feature group to one channel of one window |
| Extractor | A complete pipeline from currents to a feature list |
| Factory | Assembles all extractors for the project |

## Rules for modifying this folder

| A change here affects the database, training and inference | Rebuild the database and retrain after changing features |
| Keep feature names stable | The trained model recognises features by name |
| Prefer adding to modifying | Existing models stay valid |