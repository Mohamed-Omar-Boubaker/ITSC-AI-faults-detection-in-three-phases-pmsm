classdef BinaryHealthMetrics
    % Healthy vs faulty decision metrics (threshold on predicted severity)
    methods (Static)
        function b = compute(y, yp, thr)
            b.actualFaulty = y > 0;
            b.predictedFaulty = yp > thr;
            TP = sum(b.predictedFaulty & b.actualFaulty);
            TN = sum(~b.predictedFaulty & ~b.actualFaulty);
            FP = sum(b.predictedFaulty & ~b.actualFaulty);
            FN = sum(~b.predictedFaulty & b.actualFaulty);
            b.TP = TP; b.TN = TN; b.FP = FP; b.FN = FN;
            b.accuracy    = (TP + TN) / (TP + TN + FP + FN) * 100;
            b.precision   = TP / max(TP + FP, 1) * 100;
            b.recall      = TP / max(TP + FN, 1) * 100;
            b.specificity = TN / max(TN + FP, 1) * 100;
            b.f1 = 2 * b.precision * b.recall / max(b.precision + b.recall, eps);
            den = sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN));
            if den == 0, mcc = 0; else, mcc = (TP*TN - FP*FN) / den; end
            b.mccPct = mcc * 100;
        end
    end
end
