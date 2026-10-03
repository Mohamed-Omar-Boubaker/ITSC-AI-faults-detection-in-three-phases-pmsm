classdef FeatureFamilyCollapser
    % Collapses each feature family (e.g. I1_Peak_w1..wN) into mean (or mean+std).
    % Decisions (family size, average correlation) are fitted on TRAIN data only.
    % Name-based: transform() works on any number of windows -> usable at inference.
    properties (SetAccess = private)
        FamilyRoots = {}
        Modes = []          % 0 = keep columns as is, 1 = mean, 2 = mean + std
        TrainColumns = {}   % training column names of each family
        OutputNames = {}
    end
    properties (Constant)
        CorrThreshold = 0.9
        MinFamilySize = 3
    end
    methods
        function obj = fit(obj, X, names)
            roots = FeatureFamilyCollapser.familyRoot(names);
            [fam, ~, idx] = unique(roots, 'stable');
            modes = zeros(1, numel(fam)); cols = cell(1, numel(fam)); out = {};
            for f = 1:numel(fam)
                c = find(idx == f);
                cols{f} = names(c);
                if numel(c) < obj.MinFamilySize
                    modes(f) = 0; out = [out, names(c)];                 %#ok<AGROW>
                    continue;
                end
                R = corr(X(:, c)); n = size(R, 1);
                avgCorr = (sum(abs(R(:))) - n) / (n^2 - n);
                if avgCorr > obj.CorrThreshold
                    modes(f) = 1; out = [out, {[fam{f} '_mean']}];       %#ok<AGROW>
                else
                    modes(f) = 2; out = [out, {[fam{f} '_mean'], [fam{f} '_std']}]; %#ok<AGROW>
                end
            end
            obj.FamilyRoots = fam; obj.Modes = modes; obj.TrainColumns = cols; obj.OutputNames = out;
        end

        function Xo = transform(obj, X, names)
            roots = FeatureFamilyCollapser.familyRoot(names);
            Xo = [];
            for f = 1:numel(obj.FamilyRoots)
                if obj.Modes(f) == 0
                    [found, loc] = ismember(obj.TrainColumns{f}, names);
                    if ~all(found)
                        error('FeatureFamilyCollapser:missing', 'Missing feature(s): %s', strjoin(obj.TrainColumns{f}(~found), ', '));
                    end
                    Xo = [Xo, X(:, loc)];                                %#ok<AGROW>
                else
                    cols = find(strcmp(roots, obj.FamilyRoots{f}));
                    if isempty(cols)
                        error('FeatureFamilyCollapser:missing', 'Missing feature family: %s', obj.FamilyRoots{f});
                    end
                    Xf = X(:, cols);
                    if obj.Modes(f) == 1, Xo = [Xo, mean(Xf, 2)];                 %#ok<AGROW>
                    else,                 Xo = [Xo, mean(Xf, 2), std(Xf, 0, 2)];  %#ok<AGROW>
                    end
                end
            end
        end

        % Families whose number of windows differs from training (inference check)
        function [nDiff, example] = countMismatch(obj, names)
            roots = FeatureFamilyCollapser.familyRoot(names);
            nDiff = 0; example = '';
            for f = 1:numel(obj.FamilyRoots)
                if obj.Modes(f) == 0, continue; end
                nNow = sum(strcmp(roots, obj.FamilyRoots{f}));
                nTrain = numel(obj.TrainColumns{f});
                if nNow ~= nTrain
                    nDiff = nDiff + 1;
                    if isempty(example)
                        example = sprintf('%s: %d window(s) now vs %d in training', obj.FamilyRoots{f}, nNow, nTrain);
                    end
                end
            end
        end
    end
    methods (Static)
        function r = familyRoot(names), r = regexprep(names, '_?\d+$', ''); end
    end
end
