classdef OutOfFoldPredictor
    % Out-of-fold predictions (averaged over CV repeats) and per-fold errors
    methods (Static)
        function [oof, cnt] = compute(X, y, partitions, params, weighter)
            s = zeros(size(y)); cnt = zeros(size(y));
            for r = 1:numel(partitions)
                cvp = partitions{r};
                for k = 1:cvp.NumTestSets
                    tr = cvp.training(k); va = cvp.test(k);
                    m = SvrTrainer.train(X(tr, :), y(tr), params, weighter.compute(y(tr)));
                    s(va) = s(va) + SvrTrainer.estimate(m, X(va, :), true);
                    cnt(va) = cnt(va) + 1;
                end
            end
            oof = s ./ max(cnt, 1);
        end
        function [rmse, mae] = foldErrors(X, y, partitions, params, weighter)
            rmse = []; mae = [];
            for r = 1:numel(partitions)
                cvp = partitions{r};
                for k = 1:cvp.NumTestSets
                    tr = cvp.training(k); va = cvp.test(k);
                    m = SvrTrainer.train(X(tr, :), y(tr), params, weighter.compute(y(tr)));
                    p = SvrTrainer.estimate(m, X(va, :), true);
                    rmse(end+1) = sqrt(mean((y(va) - p).^2));            %#ok<AGROW>
                    mae(end+1)  = mean(abs(y(va) - p));                  %#ok<AGROW>
                end
            end
        end
    end
end
