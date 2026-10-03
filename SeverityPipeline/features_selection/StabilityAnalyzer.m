classdef StabilityAnalyzer
    % How often each group is picked when selection is repeated on CV training folds
    methods (Static)
        function freq = run(X, y, groups, cfg, scorer, evaluator, nonRedundant, rel, nFinal)
            SelectionRules.header('STEP 6 - FEATURE SELECTION STABILITY');
            s = rng; rng(cfg.Seed); cvp = cvpartition(size(X, 1), 'KFold', cfg.NSplits); rng(s);
            freq = containers.Map('KeyType', 'char', 'ValueType', 'double');
            for k = 1:cvp.NumTestSets
                tr = training(cvp, k); Xf = X(tr, :); yf = y(tr);
                if cfg.RecomputeRelevancePerFold
                    relF = scorer.score(Xf, yf, groups);
                    ret  = SelectionRules.relevance(relF, cfg, scorer, false);
                    pool = RedundancyRemover.remove(Xf, groups, ret, relF, cfg, false);
                else
                    pool = nonRedundant; relF = rel;
                end
                h = ForwardSelector(groups, evaluator, cfg).select(Xf, yf, pool, relF, nFinal, false);
                for i = 1:height(h)
                    g = h.group_added{i};
                    if isKey(freq, g), freq(g) = freq(g) + 1; else, freq(g) = 1; end
                end
            end
            fprintf('\nBase Feature Selection Frequency\n%s\n', repmat('-', 1, 58));
            kk = keys(freq); v = cell2mat(values(freq)); [~, o] = sort(v, 'descend');
            for i = o, fprintf('%s %d/%d = %.0f%%\n', kk{i}, v(i), cfg.NSplits, v(i)/cfg.NSplits*100); end
        end
    end
end