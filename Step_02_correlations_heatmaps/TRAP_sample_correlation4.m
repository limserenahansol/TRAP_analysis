function TRAP_sample_correlation4()
% -----------------------------------------------------------
% Build SAMPLE × SAMPLE correlation matrix using TRAP density.
%
% ALWAYS works:
%   - tolerates NaNs
%   - tolerates few samples
%   - uses bilateral averaging
%   - clusters samples (not conditions)
%
% OUTPUT:
%   dendrogram + correlation heatmap
% -----------------------------------------------------------

%% Paths
csvPath  = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";

[csvFolder,~,~] = fileparts(csvPath);
outDir = fullfile(csvFolder, "TRAP_samplecorr_output");
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf("===== SAMPLE CORRELATION ANALYSIS =====\n");

%% Load table
T = readtable(csvPath,'VariableNamingRule','preserve');

metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta = ismember(T.Properties.VariableNames, metaCols);
Node   = T(:, isMeta);

% density columns
allCols = T.Properties.VariableNames;
isDensity = contains(allCols,"density (cells/mm^3)") & ~contains(allCols,"AVERAGE");
sampleNames = allCols(isDensity);
DataTbl = T(:, sampleNames);

%% === Sample grouping ===
[GroupA, GroupB] = sampleGroups(sampleNames);

% Exclude unwanted sample (8605 black)
keep = GroupB ~= "Exclude";
sampleNames = sampleNames(keep);
GroupA = GroupA(keep);
GroupB = GroupB(keep);
DataTbl = DataTbl(:,keep);

fprintf("Using %d samples\n", numel(sampleNames));

%% === Bilateral averaging ===
acrs = string(Node.acronym);
base = erase(acrs, ["-L","-R"]);
uniq = unique(base,'stable');

Y = array2table(nan(numel(uniq), numel(sampleNames)), ...
        "VariableNames", sampleNames);
Y.region = uniq;

for u = 1:numel(uniq)
    r = uniq{u};
    idxL = acrs == r+"-L";
    idxR = acrs == r+"-R";

    if any(idxL) && any(idxR)
        Y{u,1:end-1} = (DataTbl{idxL,:} + DataTbl{idxR,:})/2;
    elseif any(idxL)
        Y{u,1:end-1} = DataTbl{idxL,:};
    elseif any(idxR)
        Y{u,1:end-1} = DataTbl{idxR,:};
    end
end

regionMat = Y{:,1:end-1};

%% === Region-wise z-scoring ===
X = zscore(regionMat,0,2);

%% === SAMPLE × SAMPLE correlation ===
C = corrcoef(X);   % Nsample × Nsample

% -------- FIX NANs (critical!) ----------
C(isnan(C)) = 0;
C(1:size(C,1)+1:end) = 1;   % diagonal = 1

%% === Convert to distance matrix ===
D = 1 - C;
D(1:size(D,1)+1:end) = 0;   % enforce valid diagonal

% -------- FIX negative distances due to rounding ----------
D(D < 0) = 0;

%% === Hierarchical clustering ===
Z = linkage(squareform(D), 'average');
[~,~,perm] = dendrogram(Z,0);

C_re = C(perm,perm);
names_re = sampleNames(perm);
GA = GroupA(perm);
GB = GroupB(perm);

%% === Figure ===
fig = figure('Color','w','Position',[200 200 900 900]);
tiledlayout(2,1,"TileSpacing","compact");

% ---- dendrogram ----
nexttile(1);
dendrogram(Z,0,'Labels',names_re,'Orientation','top');
xtickangle(45);
title("Sample dendrogram (density)");
ylabel("1 - r");

% ---- heatmap ----
nexttile(2);
imagesc(C_re, [0 1]);
axis square;
colormap(parula);
colorbar;
title("Sample × Sample correlation");
xticks(1:numel(names_re));
yticks(1:numel(names_re));
xticklabels(names_re);
yticklabels(names_re);
xtickangle(45);

% annotate Active/Passive + Phase
for i=1:numel(names_re)
    text(0.5, i, sprintf("%s | %s", GA(i), GB(i)), ...
        'HorizontalAlignment','right', ...
        'FontSize',7);
end

outfile = fullfile(outDir,"SampleCorr_density.png");
exportgraphics(fig, outfile, "Resolution",300);
close(fig);

fprintf("Saved sample correlation to:\n %s\n", outfile);
end



%% ============================================================
% CUSTOM GROUP RULES (your exact new rules)
%% ============================================================
function [GA,GB] = sampleGroups(sampleNames)
n = numel(sampleNames);
GA = strings(n,1);
GB = strings(n,1);

for i = 1:n
    nm = string(sampleNames{i});

    % Active vs Passive
    if contains(nm,"black")
        GA(i)="Passive";
    else
        GA(i)="Active";
    end

    % Phase groups (updated)
    if contains(nm,"8605_black")
        GB(i)="Exclude";
    elseif contains(nm,"7597")
        GB(i)="Withdrawal";
    elseif contains(nm,"8606_red")
        GB(i)="Reinstatement";   % moved from Reexposure to Reinst.
    elseif contains(nm,"8768_one") || contains(nm,"8606_white") || ...
           contains(nm,"8605_white") || contains(nm,"8606_black")
        GB(i)="Reinstatement";
    else
        GB(i)="Unknown";
    end
end
end
