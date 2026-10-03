projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

src = 'C:\Users\moham\OneDrive\Desktop\production\pmsm_itsc_dataset.h5';
dst = 'C:\Users\moham\OneDrive\Desktop\production\pmsm_itsc_dataset_shuffled.h5';
Hdf5Shuffler(src, dst, 42).run();