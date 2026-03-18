function trap_ensure_dir(p)
    if ~exist(p, 'dir'), mkdir(p); end
end
