function TRAP_export_depth56_region_names()
% TRAP_export_depth56_region_names
%
% TRAP_region_clusters_by_phase_density_v2 가 만든
% TRAP_downstream_input.mat 에서
%  - depth 5–6, L/R averaged 된 NodeSel 을 불러와서
%  - 1열: Allen acronym
%  - 2열: Allen full region name
% 만 뽑아서 CSV로 저장.
%
% OUTPUT:
%   <baseDir>/Depth56_region_acronym_fullname.csv
%
% baseDir 는 v2 코드에서 사용한 폴더와 동일해야 함.

    % ==== 1. v2 output 위치 설정 (필요시 경로만 수정) ====
    baseDir = "C:\Users\hsollim\Downloads\TRAP_region_clusters_by_phase_density_v2";
    matPath = fullfile(baseDir, "TRAP_downstream_input.mat");

    if ~isfile(matPath)
        error("Cannot find TRAP_downstream_input.mat at: %s\n먼저 v2 코드를 실행해서 .mat을 생성해야 합니다.", matPath);
    end

    % ==== 2. .mat 파일에서 NodeSel 불러오기 ====
    S = load(matPath);
    if ~isfield(S,'NodeSel')
        error("TRAP_downstream_input.mat 안에 NodeSel 변수가 없습니다.");
    end

    NodeSel = S.NodeSel;

    % NodeSel 은 이미:
    %  - depth 5–6 만 포함
    %  - L/R 평균 후 acronym 에서 '-L' 제거된 상태
    % 라고 가정 (v2 코드와 동일한 로직)

    acr = string(NodeSel.acronym);
    fullName = string(NodeSel.name);

    % 혹시 같은 acronym 이 depth 5, 6 두 줄로 중복되어 있을 수 있으니
    % acronym 기준으로 unique 처리 (stable 옵션)
    [~, uniqIdx] = unique(acr, 'stable');

    acr_uniq  = acr(uniqIdx);
    name_uniq = fullName(uniqIdx);

    % ==== 3. 테이블 만들어서 CSV로 저장 ====
    outTbl = table(acr_uniq, name_uniq, ...
        'VariableNames', {'Acronym','FullName'});

    outPath = fullfile(baseDir, "Depth56_region_acronym_fullname.csv");
    writetable(outTbl, outPath);

    fprintf("Exported %d unique depth 5–6 regions to:\n  %s\n", numel(acr_uniq), outPath);
end
