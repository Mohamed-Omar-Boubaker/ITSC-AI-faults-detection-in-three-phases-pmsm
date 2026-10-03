classdef ConfusionMetrics
    methods (Static)
        function ba = balancedAccuracy(C), ba = mean(diag(C) ./ max(sum(C, 2), 1)); end

        function m = multiClass(C)
            nc = size(C, 1); total = sum(C(:));
            prec = zeros(nc,1); rec = zeros(nc,1); spec = zeros(nc,1); f1 = zeros(nc,1); mcc = zeros(nc,1);
            for c = 1:nc
                tp = C(c,c); fp = sum(C(:,c)) - tp; fn = sum(C(c,:)) - tp; tn = total - tp - fp - fn;
                if (tp + fp) > 0, prec(c) = tp / (tp + fp); end
                if (tp + fn) > 0, rec(c)  = tp / (tp + fn); end
                if (tn + fp) > 0, spec(c) = tn / (tn + fp); end
                if (prec(c) + rec(c)) > 0, f1(c) = 2 * prec(c) * rec(c) / (prec(c) + rec(c)); end
                den = sqrt(double(tp+fp) * double(tp+fn) * double(tn+fp) * double(tn+fn));
                if den > 0, mcc(c) = (tp*tn - fp*fn) / den; end
            end
            n = sum(C, 2);
            m.prec = prec; m.rec = rec; m.spec = spec; m.f1 = f1; m.mcc = mcc;
            m.macro_prec = mean(prec)*100; m.macro_rec = mean(rec)*100; m.macro_spec = mean(spec)*100;
            m.macro_f1 = mean(f1)*100; m.macro_mcc = mean(mcc)*100;
            m.weight_prec = sum(prec.*n)/sum(n)*100; m.weight_rec = sum(rec.*n)/sum(n)*100;
            m.weight_spec = sum(spec.*n)/sum(n)*100; m.weight_f1 = sum(f1.*n)/sum(n)*100;
            m.weight_mcc = sum(mcc.*n)/sum(n)*100;
        end
    end
end
