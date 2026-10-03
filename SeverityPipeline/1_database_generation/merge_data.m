% File paths
dir_path = 'C:\Users\moham\OneDrive\Desktop\production';
file1_path = fullfile(dir_path, 'all_time_frequency_feat.mat');
file2_path = fullfile(dir_path, 'all_SWT_feats.mat');
output_path = fullfile(dir_path, 'Feats_DB.mat');

% Load datasets
data1 = load(file1_path);
data2 = load(file2_path);

% 1. Verify exact row-by-row alignment between files
if istable(data1.feat_table)
    n_rows1 = height(data1.feat_table);
    n_rows2 = height(data2.feat_table);
else
    n_rows1 = size(data1.feat_table, 1);
    n_rows2 = size(data2.feat_table, 1);
end

assert(n_rows1 == n_rows2, ...
    'Row mismatch error: File 1 has %d rows, but File 2 has %d rows.', n_rows1, n_rows2);

% 2. Horizontally concatenate features (feat_table) and feature names (names)
feat_table = [data1.feat_table, data2.feat_table];

if iscell(data1.names) && iscell(data2.names)
    names = [data1.names(:)', data2.names(:)'];
else
    names = [data1.names, data2.names];
end

% 3. Retain a single copy of the label vectors (degree and percent)
degree = data1.degree;
percent = data1.percent;

% 4. Save to Full_DB.mat
save(output_path, 'feat_table', 'names', 'degree', 'percent');
fprintf('Merged dataset saved successfully to: %s\n', output_path);