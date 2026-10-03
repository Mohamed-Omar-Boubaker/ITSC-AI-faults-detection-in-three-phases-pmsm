classdef FeatureVariancePlotter
    % Horizontal log-scale bar chart of feature variances (values printed beside each bar)
    methods (Static)
        function T = readTable(path)              % text file: "<feature> <variance>" per line
            fid = fopen(path); C = textscan(fid, '%s %f'); fclose(fid);
            T = table(C{1}, C{2}, 'VariableNames', {'Feature', 'Variance'});
        end
        function plot(T, outBase)
            T = sortrows(T, 'Variance', 'descend'); n = height(T);
            fprintf('Number of features: %d\n', n);
            fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 9 22]);
            ax = axes(fig); y = 1:n;
            barh(ax, y, T.Variance, 0.72, 'FaceColor', 'k', 'EdgeColor', 'k'); hold(ax, 'on');
            set(ax, 'XScale', 'log', 'YDir', 'reverse', 'YTick', y, 'YTickLabel', T.Feature, ...
                'TickLabelInterpreter', 'none', 'FontSize', 7.5, 'Box', 'off', ...
                'XGrid', 'on', 'GridLineStyle', ':', 'Layer', 'bottom');
            ylim(ax, [0.5 n + 0.5]); xlim(ax, [min(T.Variance)/2, max(T.Variance)*20]);
            text(ax, T.Variance * 1.08, y, compose('%.3e', T.Variance), 'FontSize', 5.8, ...
                'VerticalAlignment', 'middle', 'Interpreter', 'none');
            xlabel(ax, 'Pooled variance', 'FontSize', 10); ylabel(ax, 'Feature', 'FontSize', 10);
            title(ax, sprintf('Variance of the %d Features in the Initial Feature Pool', n), 'FontSize', 12);
            exportgraphics(fig, [outBase '.png'], 'Resolution', 600);
            exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
        end
    end
end