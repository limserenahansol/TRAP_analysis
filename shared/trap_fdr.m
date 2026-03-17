function q = trap_fdr(p, method)
%TRAP_FDR  BH (Benjamini–Hochberg) or BY (Benjamini–Yekutieli) adjusted q-values.

    if nargin < 2
        method = 'BH';
    end
    p = p(:);
    m = numel(p);
    q = nan(m, 1);

    if strcmpi(method, 'BY')
        cm = sum(1 ./ (1:m));
        [ps, idx] = sort(p);
        adj = nan(m, 1);
        for i = 1:m
            if isnan(ps(i))
                adj(i) = NaN;
            else
                adj(i) = min(1, ps(i) * m * cm / i);
            end
        end
        for i = m - 1:-1:1
            if ~isnan(adj(i)) && ~isnan(adj(i + 1))
                adj(i) = min(adj(i), adj(i + 1));
            end
        end
        q(idx) = adj;
        return;
    end

    % BH
    [ps, idx] = sort(p);
    rank = (1:m)';
    prev = 1;
    for i = m:-1:1
        if isnan(ps(i))
            continue;
        end
        qi = ps(i) * m / rank(i);
        if i < m
            qi = min(qi, prev);
        end
        q(idx(i)) = qi;
        prev = qi;
    end
end
