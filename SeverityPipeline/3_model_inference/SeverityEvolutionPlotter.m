classdef SeverityEvolutionPlotter
    % Top: severity per segment (+ class boundaries). Bottom: class per segment.
    methods (Static)
        function plot(res)
            K = numel(res.ClassNames);
            t = (res.TStart + res.TEnd) / 2;                         % segment centre
            figure('Name', 'Severity evolution', 'Position', [100 100 950 650], 'Color', 'w');

            subplot(2, 1, 1);
            plot(t, res.Severity, '-o', 'LineWidth', 1.8, 'MarkerFaceColor', [0.12 0.47 0.71]); hold on; grid on;
            for e = res.Edges(isfinite(res.Edges))
                yline(e, ':', 'Color', [0.5 0.5 0.5]);               % class boundaries
            end
            for k = 1:numel(t)
                text(t(k), res.Severity(k), ['  ' strtok(res.ClassNames{res.ClassIdx(k)}, ' ')], ...
                    'VerticalAlignment', 'bottom', 'FontSize', 9);
            end
            xticks(t); xlim([min(res.TStart), max(res.TEnd)]);
            ylabel('Predicted severity (%)', 'FontWeight', 'bold');
            title('Fault severity evolution (one value per segment)');

            subplot(2, 1, 2);
            stairs([res.TStart, res.TEnd(end)], [res.ClassIdx, res.ClassIdx(end)], 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); hold on;
            plot(t, res.ClassIdx, 's', 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
            grid on; yticks(1:K); yticklabels(res.ClassNames); ylim([0.5, K + 0.5]);
            xticks(t); xlim([min(res.TStart), max(res.TEnd)]);
            xlabel('Time in original signal (s)', 'FontWeight', 'bold');
            ylabel('Class', 'FontWeight', 'bold');
            title('Severity class evolution');
        end
    end
end
