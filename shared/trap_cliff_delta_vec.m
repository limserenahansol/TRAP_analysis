function d = trap_cliff_delta_vec(x, y)
%TRAP_CLIFF_DELTA_VEC  Cliff's delta for two vectors (vectorized).

    x = x(~isnan(x(:)));
    y = y(~isnan(y(:)));
    if isempty(x) || isempty(y)
        d = NaN;
        return;
    end
    x = x(:);
    y = y(:);
    gt = sum(x > y', 'all');
    lt = sum(x < y', 'all');
    d = (gt - lt) / (numel(x) * numel(y));
end
