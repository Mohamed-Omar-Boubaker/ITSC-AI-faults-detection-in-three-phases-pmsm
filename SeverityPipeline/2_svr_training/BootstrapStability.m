classdef BootstrapStability
    methods (Static)
        function r = run(y, oof, edges, nBoot)
            n = numel(y); accs = zeros(nBoot, 1);
            for b = 1:nBoot
                idx = randi(n, n, 1);
                tc = discretize(y(idx), edges); pc = discretize(oof(idx), edges);
                v = ~isnan(tc) & ~isnan(pc);
                accs(b) = mean(tc(v) == pc(v));
            end
            r.mean_acc = mean(accs)*100; r.std_acc = std(accs)*100;
            r.ci_lower = prctile(accs, 2.5)*100; r.ci_upper = prctile(accs, 97.5)*100;
        end
    end
end
