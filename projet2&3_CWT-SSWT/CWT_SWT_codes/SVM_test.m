%% =========================================================================
%  COMPLETE MATLAB PIPELINE FOR MECHANICAL FAULT SEVERITY REGRESSION
%  Publication-quality regression with SVR, Fuzzy SVR, and LS-SVR
%  using 10-fold cross-validation, PCA, and Bayesian optimization.
% =========================================================================

clear; clc; close all;
rng(42, 'twister');   % Reproducibility

fprintf('========================================================\n');
fprintf('   FAULT SEVERITY REGRESSION PIPELINE\n');
fprintf('========================================================\n\n');

%% - 1. SYNTHETIC DATA GENERATION -
% Replace this block with your actual feature matrix X and numerical
% severity vector y.
% X: [N x M] features, y: [N x 1] CONTINUOUS numerical fault severity
% (physically ordered, e.g. around 1%, 3%, 5%), NOT discrete class IDs.
N = 300;        % number of samples
M = 20;         % number of features
nLevels = 3;    % three nominal severity neighborhoods: 1%, 3%, 5%
severityLevels = [1, 3, 5];   % nominal severity values (%)

% Create feature-space centroids per nominal severity level (M features)
levelMeans = zeros(M, nLevels);
% Separate the first 3 features to distinguish severity neighborhoods
levelMeans(1:3, 1) = [0; 0; 0];
levelMeans(1:3, 2) = [3; 3; 3];
levelMeans(1:3, 3) = [6; 6; 6];
% All other features remain zero for all levels (i.e., they add only noise)

X = randn(N, M);
y = zeros(N, 1);
samplesPerLevel = floor(N / nLevels);
for c = 1:nLevels
    idx = (c-1)*samplesPerLevel + 1 : min(c*samplesPerLevel, N);
    X(idx, :) = X(idx, :) + repmat(levelMeans(:,c)', length(idx), 1);
    % Continuous numerical severity target: nominal level plus small
    % measurement-style jitter, so y is a genuine regression target
    % rather than a discrete class ID.
    y(idx) = severityLevels(c) + 0.15 * randn(length(idx), 1);
end
% Shuffle
shuf = randperm(N);
X = X(shuf, :);
y = y(shuf);

fprintf('Data dimensions: X [%d x %d], y [%d x 1]\n', size(X,1), size(X,2), length(y));
fprintf('Target severity summary: min %.3f | max %.3f | mean %.3f | std %.3f\n', ...
    min(y), max(y), mean(y), std(y));
fprintf('\n');

%% - 2. CROSS-VALIDATION SETUP -
K = 10;
% NOTE: y is now a continuous regression target, so we partition on the
% number of observations rather than stratifying by class label (the
% previous cvpartition(y,'KFold',K) call relied on y being a discrete
% grouping variable, which is no longer meaningful here).
cv = cvpartition(N, 'KFold', K);

%% - 3. PREALLOCATION FOR RESULTS -
models = {'Standard SVR', 'Fuzzy SVR', 'LS-SVR'};
nModels = 3;

% Global out-of-fold predictions and their matching true targets, per
% model (accumulated in fold-processing order so the pairing is always
% correct for the regression scatter/residual plots).
globalPred = cell(1, nModels);
globalTrue = cell(1, nModels);

% Per‑fold metrics: store each fold's metrics in a cell array to avoid
% "dissimilar structures" errors.
foldMetrics = cell(1, nModels);
for m = 1:nModels
    foldMetrics{m} = cell(K, 1);   % each cell will hold a struct
end

% Store PCA info per fold
pcaComponents = zeros(K,1);
pcaVariance   = zeros(K,1);
foldNum = 0;

% Check LS‑SVM toolbox availability (same toolbox entry point is reused
% for regression via the function-estimation type 'f').
lssvmAvailable = exist('trainlssvm', 'file') == 2;
if ~lssvmAvailable
    fprintf('LS-SVMlab toolbox not detected. LS-SVR will be skipped.\n\n');
end

%% - 4. MAIN CROSS-VALIDATION LOOP -
for k = 1:K
    foldNum = k;
    fprintf('\n=========================================\n');
    fprintf('        PROCESSING FOLD %d / %d\n', k, K);
    fprintf('=========================================\n');
    
    % Split data
    trainIdx = training(cv, k);
    testIdx  = test(cv, k);
    X_train  = X(trainIdx, :);
    y_train  = y(trainIdx);
    X_test   = X(testIdx, :);
    y_test   = y(testIdx);
    
    %  Normalization (Z-score) 
    mu    = mean(X_train);
    sigma = std(X_train);
    sigma(sigma == 0) = 1;   % avoid division by zero
    X_train_scaled = (X_train - mu) ./ sigma;
    X_test_scaled  = (X_test  - mu) ./ sigma;
    
    %  PCA on training data 
    [coeff, X_train_pca, ~, ~, explained] = pca(X_train_scaled);
    cumExplained = cumsum(explained);
    nComp = find(cumExplained >= 95, 1, 'first');
    if isempty(nComp)
        nComp = size(X_train_scaled, 2);
    end
    X_train_reduced = X_train_pca(:, 1:nComp);
    X_test_reduced  = X_test_scaled * coeff(:, 1:nComp);
    
    pcaComponents(k) = nComp;
    pcaVariance(k)   = cumExplained(nComp);
    
    fprintf('  PCA: %d components retained, %.2f%% variance explained\n', nComp, pcaVariance(k));
    
    %  Inner CV for hyperparameter optimization 
    % NOTE: y_train is continuous, so partition on sample count rather
    % than stratifying by label, as in Section 2.
    innerCv = cvpartition(numel(y_train), 'KFold', 5);
    optOptions = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                        'MaxObjectiveEvaluations', 30, ...
                        'ShowPlots', false, ...
                        'Verbose', 0, ...
                        'CVPartition', innerCv);
    
    % - Train each regressor 
    for m = 1:nModels
        % Skip LS‑SVR if not available
        if m == 3 && ~lssvmAvailable
            continue;
        end
        
        modelName = models{m};
        fprintf('\n   Training %s \n', modelName);
        
        % Measure training time
        tic;
        switch m
            case 1  % Standard SVR (RBF)
                model = fitrsvm(X_train_reduced, y_train, ...
                                 'KernelFunction', 'gaussian', ...
                                 'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale', 'Epsilon'}, ...
                                 'HyperparameterOptimizationOptions', optOptions);
                
            case 2  % Fuzzy SVR (sample weights from target-based quantile neighborhoods)
                % ADAPTATION FOR REGRESSION:
                % The original fuzzy weighting grouped samples by class
                % label, then weighted each sample by its feature-space
                % distance to its class centroid (closer to centroid =
                % more "typical" = higher weight). There are no classes
                % in regression, so instead we group training samples
                % into quantile bins of the CONTINUOUS target y_train
                % (an ordinal proxy for "local severity neighborhood"),
                % then apply the exact same centroid-distance weighting
                % logic within each bin. This preserves the original
                % fuzzy concept (down-weight feature-space outliers
                % relative to their local target neighborhood) while
                % replacing only the piece that fundamentally required
                % discrete classes.
                nFuzzyBins = 5;
                edges = quantile(y_train, linspace(0, 1, nFuzzyBins+1));
                edges(1) = -inf; edges(end) = inf;
                binIdx = discretize(y_train, edges);

                weights = zeros(size(y_train));
                for b = 1:nFuzzyBins
                    mask = (binIdx == b);
                    if sum(mask) == 0, continue; end
                    centroid = mean(X_train_reduced(mask, :), 1);
                    dist = sqrt(sum((X_train_reduced(mask, :) - centroid).^2, 2));
                    maxDist = max(dist);
                    if maxDist == 0
                        weights(mask) = 1;
                    else
                        % Weight ranges from 0.5 (farthest) to 1 (centroid)
                        weights(mask) = 1 - 0.5 * (dist / maxDist);
                    end
                end
                model = fitrsvm(X_train_reduced, y_train, ...
                                 'KernelFunction', 'gaussian', ...
                                 'Weights', weights, ...
                                 'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale', 'Epsilon'}, ...
                                 'HyperparameterOptimizationOptions', optOptions);
                
            case 3  % True LS‑SVR (requires LS‑SVMlab): function estimation ('f')
                    % instead of classification ('c'); 'mse' cost instead
                    % of 'misclass'.
                [gam, sig2] = tunelssvm({X_train_reduced, y_train, 'f', [], [], 'RBF_kernel'}, ...
                                        'simplex', 'crossvalidatelssvm', {10, 'mse'});
                model = trainlssvm({X_train_reduced, y_train, 'f', gam, sig2, 'RBF_kernel'});
        end
        trainTime = toc;
        
        % Prediction
        tic;
        switch m
            case {1, 2}
                y_pred = predict(model, X_test_reduced);
            case 3
                y_pred = simlssvm(model, X_test_reduced);
        end
        predTime = toc;
        
        y_pred = double(y_pred);   % ensure numeric, continuous regression output
                                    % (never rounded/thresholded to a class label)
        
        % Compute regression metrics using helper function
        metrics = computeRegressionMetrics(y_test, y_pred);
        metrics.trainTime = trainTime;
        metrics.predTime  = predTime;
        metrics.nComp = nComp;
        metrics.varExplained = pcaVariance(k);
        
        % Store out-of-fold predictions AND their matching true targets
        globalPred{m} = [globalPred{m}; y_pred];
        globalTrue{m} = [globalTrue{m}; y_test];
        % Store per‑fold metrics in cell array
        foldMetrics{m}{k} = metrics;
        
        % Print per‑fold summary (key metrics)
        fprintf('    RMSE: %.4f\n', metrics.rmse);
        fprintf('    MAE: %.4f\n', metrics.mae);
        fprintf('    R^2: %.4f\n', metrics.r2);
        fprintf('    Correlation: %.4f\n', metrics.corr);
        fprintf('    Training time: %.4f s\n', metrics.trainTime);
        fprintf('    Prediction time: %.4f s\n', metrics.predTime);
    end
end

%% - 5. AGGREGATE RESULTS -
fprintf('\n\n======================================================\n');
fprintf('     FINAL AGGREGATED PERFORMANCE REPORT\n');
fprintf('======================================================\n');

% Define metrics to summarize
metricNames = {'rmse', 'mse', 'mae', 'r2', 'corr', 'nrmse', 'trainTime', 'predTime'};

meanStd = cell(1, nModels);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        meanStd{m} = [];
        continue;
    end
    nFolds = K;
    vals = zeros(nFolds, length(metricNames));
    for k = 1:nFolds
        met = foldMetrics{m}{k};   % extract struct from cell
        for i = 1:length(metricNames)
            if isfield(met, metricNames{i})
                vals(k, i) = met.(metricNames{i});
            else
                vals(k, i) = NaN;
            end
        end
    end
    meanVals = nanmean(vals, 1);
    stdVals  = nanstd(vals, 0, 1);
    meanStd{m} = table(meanVals', stdVals', ...
                       'VariableNames', {'Mean', 'Std'}, ...
                       'RowNames', metricNames);
end

% Print comparison table (mean ± std) for key metrics
keyMetrics = {'rmse', 'mae', 'r2', 'corr', 'trainTime', 'predTime'};
fprintf('\nComparison Table (Mean ± Std):\n');
fprintf('%-20s', 'Metric');
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        fprintf('  %-20s', 'LS-SVR (skip)');
    else
        fprintf('  %-20s', models{m});
    end
end
fprintf('\n');
for i = 1:length(keyMetrics)
    fprintf('%-20s', keyMetrics{i});
    for m = 1:nModels
        if m == 3 && ~lssvmAvailable
            fprintf('  %-20s', 'N/A');
        else
            tbl = meanStd{m};
            mu = tbl.Mean(keyMetrics{i});
            sd = tbl.Std(keyMetrics{i});
            fprintf('  %6.4f ± %6.4f', mu, sd);
        end
    end
    fprintf('\n');
end

%% - 6. GLOBAL OUT-OF-FOLD PREDICTIONS -
% globalPred{m} / globalTrue{m} already hold the complete set of
% out-of-fold predictions and their matching true severities for each
% model (accumulated fold-by-fold in Section 4), ready for the
% regression diagnostic plots below.

%% - 7. PLOTS -
% Figure 1: PCA cumulative explained variance (from last fold)
figure('Name', 'PCA Explained Variance', 'Position', [100 100 800 600]);
plot(cumExplained, 'o-', 'LineWidth', 2);
xlabel('Number of Principal Components');
ylabel('Cumulative Explained Variance (%)');
title('PCA Cumulative Explained Variance');
grid on;
ylim([0 100]);
hold on;
line([nComp nComp], ylim, 'Color', 'r', 'LineStyle', '');
text(nComp+0.5, 50, sprintf('95%% at %d comps', nComp), 'Color', 'r');
hold off;

% Figure 2: Predicted vs. true severity (scatter + y=x reference) and residuals
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        continue;
    end
    yt = globalTrue{m};
    yp = globalPred{m};
    resid = yp - yt;

    figure('Name', ['Regression Fit - ' models{m}], 'Position', [100 100 1200 500]);

    subplot(1,2,1);
    scatter(yt, yp, 25, 'filled'); hold on;
    lims = [min([yt; yp]), max([yt; yp])];
    plot(lims, lims, 'r--', 'LineWidth', 1.5);
    xlabel('True Severity'); ylabel('Predicted Severity');
    title(['Predicted vs. True - ' models{m}]);
    grid on; axis equal; xlim(lims); ylim(lims);
    hold off;

    subplot(1,2,2);
    scatter(yt, resid, 25, 'filled'); hold on;
    yline(0, 'r--', 'LineWidth', 1.5);
    xlabel('True Severity'); ylabel('Residual (Predicted - True)');
    title(['Residuals - ' models{m}]);
    grid on;
    hold off;
end

% Figure 3: Residual distribution per model (replaces ROC curves)
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        continue;
    end
    resid = globalPred{m} - globalTrue{m};

    figure('Name', ['Residual Distribution - ' models{m}], 'Position', [100 100 800 600]);
    histogram(resid, 20, 'Normalization', 'pdf'); hold on;
    xlabel('Residual (Predicted - True)'); ylabel('Density');
    title(sprintf('Residual Distribution - %s\nMean = %.4f, Std = %.4f', ...
                  models{m}, mean(resid), std(resid)));
    grid on;
    hold off;
end

% Figure 4: Bar chart comparing regressors (selected metrics)
metricsToPlot = {'rmse', 'mae', 'r2', 'corr'};
nMetrics = length(metricsToPlot);
means = zeros(nModels, nMetrics);
stds = zeros(nModels, nMetrics);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        means(m,:) = NaN;
        stds(m,:) = NaN;
        continue;
    end
    tbl = meanStd{m};
    for i = 1:nMetrics
        means(m,i) = tbl.Mean(metricsToPlot{i});
        stds(m,i)  = tbl.Std(metricsToPlot{i});
    end
end

figure('Name', 'Performance Comparison', 'Position', [100 100 1000 600]);
bar(means');
hold on;
% Add error bars
for i = 1:nMetrics
    x = (1:nModels) + (i-1)*0.2 - 0.3;   % approximate positions
    errorbar(x, means(:,i), stds(:,i), 'k.', 'LineWidth', 1);
end
hold off;
set(gca, 'XTick', 1:nModels, 'XTickLabel', models);
xlabel('Regressor');
ylabel('Metric Value');
title('Comparison of Regressors');
legend(metricsToPlot, 'Location', 'best');
grid on;

% Figure 5: Boxplots of fold RMSE
foldRMSE = zeros(K, nModels);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        foldRMSE(:,m) = NaN;
        continue;
    end
    for k = 1:K
        foldRMSE(k,m) = foldMetrics{m}{k}.rmse;
    end
end
figure('Name', 'Fold RMSE', 'Position', [100 100 800 600]);
boxplot(foldRMSE, 'Labels', models);
ylabel('RMSE');
title('Distribution of Fold RMSE');
grid on;

% Figure 6: Radar chart (spider plot) comparing regressors
% NOTE: RMSE/MAE are lower-is-better while R2/correlation are
% higher-is-better. To keep a single consistent "higher = better" radar
% axis, error metrics are inverted before normalization.
radarMetrics = {'rmse', 'mae', 'r2', 'corr'};
lowerIsBetter = [true, true, false, false];
nRadar = length(radarMetrics);
radarData = zeros(nModels, nRadar);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        radarData(m,:) = NaN;
        continue;
    end
    tbl = meanStd{m};
    for i = 1:nRadar
        radarData(m,i) = tbl.Mean(radarMetrics{i});
    end
end

normRadar = zeros(size(radarData));
for i = 1:nRadar
    col = radarData(:,i);
    if lowerIsBetter(i)
        col = -col;   % invert so higher = better, matching the other metrics
    end
    colMax = max(col);
    colMin = min(col);
    if colMax == colMin
        normRadar(:,i) = 1;
    else
        normRadar(:,i) = (col - colMin) / (colMax - colMin);
    end
end

theta = linspace(0, 2*pi, nRadar+1);
theta = theta(1:end-1);
figure('Name', 'Radar Chart', 'Position', [100 100 800 600]);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        continue;
    end
    r = [normRadar(m,:), normRadar(m,1)];
    polarplot([theta, theta(1)], r, 'o-', 'LineWidth', 2);
    hold on;
end
hold off;
ax = gca;
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = radarMetrics;
legend(models(~strcmp(models,'LS-SVR') | lssvmAvailable), 'Location', 'best');
title('Radar Chart of Regressor Performance (normalized, higher = better)');

% Figure 7: Training time comparison
trainTimes = zeros(K, nModels);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        trainTimes(:,m) = NaN;
        continue;
    end
    for k = 1:K
        trainTimes(k,m) = foldMetrics{m}{k}.trainTime;
    end
end
figure('Name', 'Training Times', 'Position', [100 100 800 600]);
bar(mean(trainTimes,1));
hold on;
errorbar(1:nModels, mean(trainTimes,1), std(trainTimes,1), 'k.', 'LineWidth', 1);
hold off;
set(gca, 'XTick', 1:nModels, 'XTickLabel', models);
ylabel('Training Time (s)');
title('Average Training Time per Regressor');
grid on;

% Figure 8: Prediction time comparison
predTimes = zeros(K, nModels);
for m = 1:nModels
    if m == 3 && ~lssvmAvailable
        predTimes(:,m) = NaN;
        continue;
    end
    for k = 1:K
        predTimes(k,m) = foldMetrics{m}{k}.predTime;
    end
end
figure('Name', 'Prediction Times', 'Position', [100 100 800 600]);
bar(mean(predTimes,1));
hold on;
errorbar(1:nModels, mean(predTimes,1), std(predTimes,1), 'k.', 'LineWidth', 1);
hold off;
set(gca, 'XTick', 1:nModels, 'XTickLabel', models);
ylabel('Prediction Time (s)');
title('Average Prediction Time per Regressor');
grid on;

fprintf('\nAll plots generated.\n');
fprintf('Pipeline complete.\n');

%% ======================== HELPER FUNCTIONS ==============================

function metrics = computeRegressionMetrics(y_true, y_pred)
    % Regression counterpart of the former computeAllMetrics: returns
    % RMSE, MSE, MAE, R^2, Pearson correlation, and normalized RMSE
    % between true and predicted continuous severity values.
    % trainTime / predTime / nComp / varExplained are attached by the
    % caller in the main loop, exactly as in the original pipeline
    % (computeAllMetrics also did not compute those fields itself).
    y_true = double(y_true(:));
    y_pred = double(y_pred(:));

    err = y_true - y_pred;

    mse  = mean(err.^2);
    rmse = sqrt(mse);
    mae  = mean(abs(err));

    ssRes = sum(err.^2);
    ssTot = sum((y_true - mean(y_true)).^2);
    if ssTot == 0
        r2 = NaN;
    else
        r2 = 1 - ssRes / ssTot;
    end

    if std(y_true) == 0 || std(y_pred) == 0
        corrVal = NaN;
    else
        R = corrcoef(y_true, y_pred);
        corrVal = R(1,2);
    end

    % Normalized RMSE (optional), normalized by the range of the true target
    rangeY = max(y_true) - min(y_true);
    if rangeY == 0
        nrmse = NaN;
    else
        nrmse = rmse / rangeY;
    end

    metrics = struct();
    metrics.rmse  = rmse;
    metrics.mse   = mse;
    metrics.mae   = mae;
    metrics.r2    = r2;
    metrics.corr  = corrVal;
    metrics.nrmse = nrmse;
end