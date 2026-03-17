function trap_write_folder_readme(figDir, titleStr, bodyText)
%TRAP_WRITE_FOLDER_README  00_README_THIS_FOLDER.txt in a figures directory.

    p = fullfile(figDir, '00_README_THIS_FOLDER.txt');
    fid = fopen(p, 'w');
    if fid > 0
        fprintf(fid, '%s\n\n%s\n', titleStr, bodyText);
        fclose(fid);
    end
end
