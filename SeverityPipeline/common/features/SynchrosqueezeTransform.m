classdef SynchrosqueezeTransform
    % Frequency reassignment (SWT) and bilateral time-frequency map assembly
    methods (Static)
        function Tx = reassign(W, dW, f, scales, threshRatio, da)
            f = f(:); scales = scales(:);
            [Ns, N] = size(W);
            [wSorted, sortIdx] = sort(2*pi*f, 'ascend');
            rowLookup = (1:Ns)';
            thresh = threshRatio * max(abs(W), [], 'all');
            valid = abs(W) > thresh;
            validIdx = find(valid);
            if isempty(validIdx), Tx = zeros(Ns, N, 'like', W); return; end
            IFrad = zeros(Ns, N);
            IFrad(valid) = abs(imag(dW(valid) ./ W(valid)));
            rows = interp1(wSorted, rowLookup(sortIdx), IFrad(validIdx), 'nearest', 'extrap');
            rows = max(1, min(Ns, rows));
            scalesMat = repmat(scales, 1, N); daMat = repmat(da, 1, N);
            weights = W(validIdx) .* (scalesMat(validIdx) .^ (-1)) .* daMat(validIdx);
            [~, cols] = ind2sub([Ns, N], validIdx);
            Tx = accumarray([rows, cols], weights, [Ns, N]);
        end

        function [TxAll, fAxis, TxAllC] = buildBilateral(TxPos, TxNeg, f, Cpsi)
            NumTime = size(TxPos, 2);
            fAxis = [-f; 0; flipud(f)];
            [fAxis, sortIdx] = sort(fAxis, 'ascend');
            TxAll  = [abs(TxNeg); zeros(1, NumTime); flipud(abs(TxPos))] / Cpsi;
            TxAll  = TxAll(sortIdx, :);
            TxAllC = [TxNeg; zeros(1, NumTime); flipud(TxPos)];
            TxAllC = TxAllC(sortIdx, :);
        end
    end
end
