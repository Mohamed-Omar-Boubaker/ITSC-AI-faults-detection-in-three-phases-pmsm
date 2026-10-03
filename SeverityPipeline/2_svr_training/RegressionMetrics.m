classdef RegressionMetrics
    methods (Static)
        function m = compute(y, yp)
            err = y - yp; m.err = err;
            m.rmse = sqrt(mean(err.^2));
            m.mae = mean(abs(err));
            m.medae = median(abs(err));
            m.maxerr = max(abs(err));
            m.mapeMask = abs(y) > 1e-6;
            m.mape = mean(abs(err(m.mapeMask) ./ y(m.mapeMask))) * 100;
            ratio = err ./ y; ratio(isnan(ratio)) = 0;
            m.maape = mean(atan(abs(ratio))) * (2/pi) * 100;
            m.r2 = 1 - sum(err.^2) / sum((y - mean(y)).^2);
        end
        function s = severityBreakdown(y, yp, thr)
            err = y - yp;
            h = y == 0; f = y > 0; lo = y > 0 & y <= thr; hi = y > thr;
            s.healthy_rmse = sqrt(mean(err(h).^2));  s.healthy_n = sum(h);
            s.faulty_rmse  = sqrt(mean(err(f).^2));  s.faulty_n  = sum(f);
            s.low_rmse     = sqrt(mean(err(lo).^2)); s.low_n     = sum(lo);
            s.high_rmse    = sqrt(mean(err(hi).^2)); s.high_n    = sum(hi);
        end
    end
end
