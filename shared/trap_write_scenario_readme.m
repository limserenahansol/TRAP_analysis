function trap_write_scenario_readme(dirRoot, figDir, body)
    p = fullfile(dirRoot, 'README_Scenario.txt');
    fid = fopen(p, 'w');
    if fid > 0
        fprintf(fid, '%s\n\nTables: tables\\\nFigures + .txt: figures_described\\\n', body);
        fclose(fid);
    end
    trap_write_folder_readme(figDir, 'Phase-specific Active vs Passive', ...
        [body char(10) char(10) 'Each PNG has matching .txt. CSVs in tables\.']);
end
