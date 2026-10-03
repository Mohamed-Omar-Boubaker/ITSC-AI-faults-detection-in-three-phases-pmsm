classdef BoundaryOptimizerDP
    % Optimal 1-D segmentation of the fault region by dynamic programming:
    % maximise diagonal hits + soft width bonus under min/max class-size constraints.
    properties (SetAccess = immutable)
        Eval
        TargetWidth
        WidthBonus
        NGrid
    end
    methods
        function obj = BoundaryOptimizerDP(ev, targetWidth, widthBonus, nGrid)
            obj.Eval = ev; obj.TargetWidth = targetWidth; obj.WidthBonus = widthBonus; obj.NGrid = nGrid;
        end

        function [best, found] = optimize(obj, yTrain, oof, K, minSize, maxSize, upperLimit)
            hb = obj.Eval.HealthyBound; found = false;
            best = struct('boundaries', [], 'edges', [], 'accuracy', NaN, 'balanced_accuracy', NaN, ...
                'macro_f1', NaN, 'macro_mcc', NaN, 'weight_f1', NaN, 'counts', [], 'score', -Inf);
            m = K - 2;
            faultMask = yTrain > hb; nFault = sum(faultMask);
            if nFault < (m+1) * 3, return; end
            yf = yTrain(faultMask); pf = oof(faultMask);
            maxVal = max([upperLimit; yf(:)]);
            if ~(maxVal > hb), return; end

            if m == 0
                res = obj.Eval.evaluate(yTrain, oof, [], K, minSize, maxSize);
                if ~res.feasible, return; end
                best = BoundaryOptimizerDP.pack([], res); found = true; return;
            end

            raw = linspace(hb, maxVal, obj.NGrid + 2);
            cand = unique(raw(2:end-1)); L = numel(cand);
            if L < m, return; end
            edgesAll = [hb, cand, Inf];
            binTrue = discretize(yf, edgesAll); binPred = discretize(pf, edgesAll);
            valid = ~isnan(binPred);
            trueCount = accumarray(binTrue(:), 1, [L+1, 1]);
            joint = accumarray([binTrue(valid), binPred(valid)], 1, [L+1, L+1]);
            cumTrue = [0; cumsum(trueCount)];
            intH = zeros(L+2, L+2);
            intH(2:end, 2:end) = cumsum(cumsum(joint, 1), 2);          % 2-D prefix sums
            edgeVal = [hb, cand];
            blockTP    = @(p, q) intH(q+1,q+1) - intH(p+1,q+1) - intH(q+1,p+1) + intH(p+1,p+1);
            blockCount = @(p, q) cumTrue(q+1) - cumTrue(p+1);

            NEG = -1e18; dp = NEG * ones(m+1, L+1); par = zeros(m+1, L+1); dp(1,1) = 0;
            for j = 1:m
                for kp = j:L
                    bestVal = NEG; bestK = -1;
                    for k = max(j-1, 0):(kp-1)
                        if dp(j, k+1) <= NEG/2, continue; end
                        cnt = blockCount(k, kp);
                        if cnt < minSize || cnt > maxSize, continue; end
                        width = edgeVal(kp+1) - edgeVal(k+1);
                        bonus = obj.WidthBonus * nFault * min(width, obj.TargetWidth) / obj.TargetWidth;
                        val = dp(j, k+1) + blockTP(k, kp) + bonus;
                        if val > bestVal, bestVal = val; bestK = k; end
                    end
                    dp(j+1, kp+1) = bestVal; par(j+1, kp+1) = bestK;
                end
            end

            bestTotal = NEG; bestKFinal = -1;
            for k = m:L                                                 % open-ended last class, no bonus
                if dp(m+1, k+1) <= NEG/2, continue; end
                cnt = blockCount(k, L+1);
                if cnt < minSize || cnt > maxSize, continue; end
                val = dp(m+1, k+1) + blockTP(k, L+1);
                if val > bestTotal, bestTotal = val; bestKFinal = k; end
            end
            if bestKFinal < 0, return; end

            bnd = zeros(1, m); k = bestKFinal;
            for j = m:-1:1
                bnd(j) = edgeVal(k+1); k = par(j+1, k+1);
            end
            res = obj.Eval.evaluate(yTrain, oof, bnd, K, minSize, maxSize);
            if ~res.feasible, return; end
            best = BoundaryOptimizerDP.pack(bnd, res); found = true;
        end
    end
    methods (Static)
        function b = pack(bnd, res)
            b.boundaries = bnd; b.edges = res.edges; b.accuracy = res.accuracy;
            b.balanced_accuracy = res.balanced_accuracy; b.macro_f1 = res.macro_f1;
            b.macro_mcc = res.macro_mcc; b.weight_f1 = res.weight_f1;
            b.counts = res.counts; b.score = res.score;
        end
    end
end
