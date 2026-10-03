%% ========================================================================
%  RUN THIS FILE TO BUILD THE FEATURE DATABASE (HDF5 -> New_DATABASE.mat)
%  Every feature is extracted cyclically:
%    - pipeline 1: per 0.2s window (I1/I2/CUF/CSV time, freq, geometric, nonlinear)
%    - pipeline 2: per 0.2 s window (SWT rec. energy, wavelet entropy, sub-band ratio)
%% ========================================================================
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

%% ---- User settings ----
h5file    = 'C:\Users\moham\OneDrive\Desktop\ITSC AI faults detection\SeverityPipeline\1_database_generation\pmsm_itsc_dataset_shuffled.h5';
outfile   = 'C:\Users\moham\OneDrive\Desktop\ITSC AI faults detection\SeverityPipeline\1_database_generation\FEATURES_DATABASE.mat';
fsOrig    = 100000;    % sampling frequency stored in the HDF5
fsTarget  = 10000;     % both pipelines run at this rate
startIdx  = 1;         % >1 resumes an interrupted run (reloads the existing outfile)

%% ---- Assemble and run ----
source     = Hdf5SampleSource(h5file, fsOrig);
extractors = FeatureExtractorFactory.createDefault(fsTarget);   % which features -> one place
writer     = MatDatabaseWriter(outfile);

FeatureDatabaseBuilder(source, extractors, writer, fsTarget).run(startIdx);
