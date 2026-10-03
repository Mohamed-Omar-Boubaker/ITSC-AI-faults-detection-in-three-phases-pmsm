classdef ConfusionPlotter
    methods (Static)
        % Ordinal categoricals keep the class order; plot + raw-count heatmap + console print
        function [C, metrics] = suite(trueIdx, predIdx, labels, titleStr)
            K = numel(labels);
            tc = categorical(trueIdx, 1:K, labels, 'Ordinal', true);
            pc = categorical(predIdx, 1:K, labels, 'Ordinal', true);
            C = confusionmat(trueIdx, predIdx, 'Order', 1:K);
            metrics = ConfusionMetrics.multiClass(C);

            figure('Name', [titleStr ' - Plot']);
            confusionchart(tc, pc, 'Title', [titleStr ' (Plot)'], ...
                'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
            figure('Name', [titleStr ' - Raw Values']);
            h = heatmap(labels, labels, C, 'Title', [titleStr ' (Raw Counts)']);
            h.XLabel = 'Predicted'; h.YLabel = 'True';

            fprintf('\n--- %s : Confusion Matrix (rows=true, cols=predicted) ---\n', titleStr);
            fprintf('%15s', ''); fprintf('%15s', labels{:}); fprintf('\n');
            for r = 1:K
                fprintf('%15s', labels{r}); fprintf('%15d', C(r, :)); fprintf('\n');
            end
        end
    end
end
