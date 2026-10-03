# SeverityPipeline — Estimating Motor Fault Severity from Currents

## 1. Why this project?

Electric motors (here, three-phase permanent magnet synchronous motors, PMSM) wear out. One of the most common electrical problems is a **short circuit inside the stator winding**: a few turns of copper touch each other and a small unwanted current loop appears.

At the beginning, the fault is small and the motor still works. If nobody notices, it grows, heats the motor and can destroy it. Detecting "fault / no fault" is useful, but **knowing how severe the fault is** allows maintenance to be planned before a breakdown.

This project answers one question:

> *Given the three phase currents of a motor, how severe is the stator fault, and in which severity class does it fall?*

It uses only the currents, which are cheap to measure and need no sensor on the motor itself.

## 2. The idea in one picture

```
 Three phase currents      Numbers that describe        A learned model          Severity (%)
 (ia, ib, ic)         ──►  the signal's behaviour  ──►  (Support Vector     ──►  + severity class
                           ("features")                  Regression)
```

A healthy motor draws three currents that are almost perfectly balanced. A winding fault breaks this balance and slightly distorts the waveforms. These changes are tiny and hard to see by eye, so we:

1. **Measure** the signal in many different ways (features).
2. **Teach** a machine-learning model the link between those features and the known severity.
3. **Use** the trained model on new signals.

## 3. Methodology in simple words

| Concept | Plain explanation |
|---|---|
| **Symmetrical components** | Any unbalanced three-phase system can be split into a "perfect" part (positive sequence) and an "imbalance" part (negative sequence). A fault makes the imbalance part grow. |
| **Space vector** | The three currents are combined into one rotating arrow (a complex number). Its shape reveals the fault. |
| **Cyclic analysis** | The signal is cut into consecutive short windows and every feature is computed on each window. This follows how the signal evolves instead of giving one global number. |
| **Time-domain features** | Statistics of the waveform: peak, spread, shape, "spikiness"… |
| **Frequency-domain features** | Where the energy sits in the spectrum, and how much is in harmonics. |
| **Geometric features** | Shape of the orbit drawn by the space vector (an ellipse when the system is unbalanced). |
| **Nonlinear features** | How complex or irregular the signal is (entropy, fractal dimension…). |
| **Wavelet / synchrosqueezing features** | A time–frequency "microscope" that follows the energy around the frequency where the fault leaves its signature. |
| **SVR (Support Vector Regression)** | A model that predicts a continuous value (the severity %) by comparing a new signal with the most informative training examples ("support vectors"). |
| **Severity classes** | The continuous severity is converted to a few ordered classes (healthy, mild, … severe) using boundaries found during training. |

## 4. Global architecture

```
  HDF5 raw signals
        │
        ▼
┌────────────────────────┐    feature table     ┌──────────────────────────┐   model bundle   ┌───────────────────────────┐
│ 1_database_generation  │ ───────────────────► │ 2_svr_training           │ ───────────────► │ 3_model_inference         │
│ currents → features    │   (.mat database)    │ features → trained model │   (.mat file)    │ new signal → severity     │
└───────────┬────────────┘                      └────────────┬─────────────┘                  └─────────────┬─────────────┘
            │                                                │                                              │
            └────────────────────────────┬───────────────────┴──────────────────────────────────────────────┘
                                         ▼
                              common/  (shared building blocks:
                              signal tools, feature extractors,
                              preprocessing, SVR wrapper)
```

The same feature extractors are used to build the database and to analyse a new signal. This is essential: the model only understands features computed exactly the way it was trained.

## 5. Content of the project

```
SeverityPipeline/
├── README.md                    ← this file
├── common/                      ← code shared by several parts
├── 1_database_generation/       ← build the feature database
├── 2_svr_training/              ← train and save the model
└── 3_model_inference/           ← use the model on a new signal
```

| Folder | Role | File to run |
|---|---|---|
| `1_database_generation` | Reads raw currents, extracts features cyclically, stores a table | `main_buildDatabase.m` |
| `2_svr_training` | Prepares the data, trains the SVR, defines classes, evaluates, saves a bundle | `main_trainSvr.m` |
| `3_model_inference` | Loads the bundle and a signal, splits it into segments, outputs severity and class per segment, plots the evolution | `main_classifySignal.m` |
| `common` | Shared classes (no file to run) | — |

Each folder has its own `README.md` with the details.

## 6. Dataset

The study is based on a public dataset of three-phase PMSMs with seeded stator faults:



| **Title** | *Vibration and Current Dataset of Three-Phase Permanent Magnet Synchronous Motors with Stator Faults* |
| **Authors** | W. Jung, S.-H. Yun, Y.-S. Lim, S. Cheong, Y.-H. Park (KAIST, South Korea) |
| **Article** | Data in Brief 47 (2023) 108952 — https://doi.org/10.1016/j.dib.2023.108952 |
| **Data (Mendeley Data)** | https://data.mendeley.com/datasets/rgn5brrgrn/5 |
| **DOI** | 10.17632/rgn5brrgrn.5 |
| **Licence** | CC BY 4.0 |

What it contains, in short:

- Several motors of different power, tested under the same speed and load.
- Two kinds of seeded stator faults: **inter-turn** and **inter-coil** short circuits, each at several severity levels, plus a healthy reference.
- Severity is a *fault ratio*: the share of current flowing through an artificial bypass resistor.
- Vibration and three-phase current recordings (stored as TDMS files).

This project uses **only the three-phase currents**, read from a prepared HDF5 file.

## 7. How to use it

| 1 | Run `1_database_generation/main_buildDatabase.m` to create the feature database. |
| 2 | (Optional) select/filter features in the database, keeping the feature names. |
| 3 | Run `2_svr_training/main_trainSvr.m` to train and save the model bundle. |
| 4 | Run `3_model_inference/main_classifySignal.m` on any signal stored in a `.mat` file 

=========================================================================================================================================
(Note: still nt complete as DATASET full sescription is not complete to generate sythetic samples otherwise the model will act as a KNN and predict always the wrong class "SORRY XD")
=========================================================================================================================================


Paths and settings are edited at the top of each `main_*.m` file.

## 8. Design principles (SOLID, explained simply)

| Principle | In this project |

| Single responsibility | Each class does one job (read, resample, extract, train, plot…). |
| Open for extension | A new feature is one new line in a list; a new pipeline is one new class. |
| Substitutable parts | Data sources, writers and extractors share small common interfaces. |
| Small interfaces | A class only needs to implement what it really uses. |
| Depend on abstractions | The orchestrators work with interfaces, not with specific file formats. |

## 9. Limitations to keep in mind

- The model learns from the dataset it was given: signals very different from it (other motor type, other acquisition chain) may not be handled correctly.
- Features are computed in a way tied to the training data; the signal must be prepared the same way.
- The severity is an estimate; the class is only as reliable as that estimate.

## 10. Glossary

| Term | Meaning |

| PMSM | Permanent magnet synchronous motor |
| Stator | Fixed part of the motor, carrying the windings |
| Inter-turn / inter-coil fault | Short circuit between neighbouring turns / coils |
| Feature | A number summarising a property of a signal |
| SVR | Support Vector Regression |
| Support vector | A training example that defines the model's decision |
| Cross-validation | Re-training on different data splits to test reliability |
| Bundle | The saved file containing the model and all its preprocessing steps |