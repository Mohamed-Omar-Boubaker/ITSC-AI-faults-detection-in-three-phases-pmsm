classdef GeometricFeatures < FeatureGroup
    % Orbit (ellipse) features of the complex space vector
    methods
        function r = registry(~)
            r = {'Eccentricity',    @GeometricFeatures.eccentricity; ...
                 'Orbit_Thickness', @GeometricFeatures.orbitThickness; ...
                 'Enclosed_Area',   @GeometricFeatures.enclosedArea};
        end
    end
    methods (Static)
        function [a, b] = ellipseAxes(z)
            Z = [real(z(:)), imag(z(:))];
            [~, S, ~] = svd(Z - mean(Z, 1), 0);
            sg = diag(S) * sqrt(2 / length(z));
            a = max(sg); b = min(sg);
        end
        function v = eccentricity(z), [a, b] = GeometricFeatures.ellipseAxes(z); v = sqrt(1 - (b / (a + eps))^2); end
        function v = orbitThickness(z), r = abs(z(:)); v = std(r, 1) / (mean(r) + eps); end
        function v = enclosedArea(z), [a, b] = GeometricFeatures.ellipseAxes(z); v = pi * a * b; end
    end
end
