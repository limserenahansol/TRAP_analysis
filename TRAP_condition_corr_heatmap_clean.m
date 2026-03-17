function TRAP_condition_corr_heatmap_clean()

%% === Paths ===
csvPath = "C:\Users\hsollim\Downloads\Hansol Lim density channel 561_all.csv";
T = readtable(csvPath,'VariableNamingRule','preserve');

%% === Extract meta + density columns ===
metaCols = {'id','name','acronym','parent_structure_id','depth'};
isMeta   = ismember(T.Properties.VariableNames, metaCols);
NodeFull = T(:, isMeta);

allVars = T.Properties.VariableNames;
isDens  = contains(allVars,"density (cells/mm^3)") & ~contains(allVars,"AVERAGE");
sampleNames = string(allVars(isDens));
Data = T{:, isDens};   % regions × samples

%% === Build Group labels (your final desired 4 groups) ===
n = numel(sampleNames);
Group = strings(n,1);

for i = 1:n
    nm = sampleNames(i);

    % Delivery (Active / Passive)
    isPassive = contains(nm,"black");
    delivery  = "Active";
    if isPassive, delivery = "Passive"; end

    % Phase (Withdrawal / Reinstatement)
    if contains(nm,"7597")
        phase = "Withdrawal";
    else
        phase = "Reinstatement";
    end

    % Combine
    Group(i) = delivery + "-" + phase;
end

% Exclude sample 8605_black (your rule)
badIdx = contains(sampleNames,"8605_black");
Group(badIdx) = "Exclude";

keep = Group ~= "Exclude";
sampleNames = sampleNames(keep);
Group       = Group(keep);
Data        = Data(:,keep);

%% === Bilateral averaging ===
ac = string(NodeFull.acronym);
isLeft  = endsWith(ac,"-L");
isRight = endsWith(ac,"-R");
isGlobal = ~(isLeft | isRight);

keepMask = isLeft | isGlobal;
Node = NodeFull(keepMask,:);
idxKeep = find(keepMask);
nReg = height(Node);

densMean = nan(nReg, numel(sampleNames));

for ii = 1:nReg
    idxL = idxKeep(ii);
    nmL = ac(idxL);

    if endsWith(nmL,"-L")
        base = extractBefore(nmL,"-L");
        nmR  = base+"-R";
        idxR = find(ac == nmR,1);

        if ~isempty(idxR)
            densMean(ii,:) = (Data(idxL,:) + Data(idxR,:))/2;
        else
            densMean(ii,:) = Data(idxL,:);
        end
    else
        densMean(ii,:) = Data(idxL,:);
    end
end

%% === 4-condition mean vectors ===
conds = unique(Group,'stable');
nCond = numel(conds);

CondMat = nan(nReg, nCond);

for c = 1:nCond
    idx = Group == conds(c);
    CondMat(:,c) = mean(densMean(:,idx),2,'omitnan');
end

%% === Correlation matrix (conditions × conditions) ===
C = corrcoef(CondMat);

%% === Heatmap ===
figure('Color','w','Position',[300 200 900 800]);
imagesc(C,[0 1]);
axis square;
colormap(jet);
colorbar;

xticks(1:nCond); yticks(1:nCond);
xticklabels(conds); yticklabels(conds);
xtickangle(45);

title('Condition correlation (Delivery × Phase, density)','FontSize',16);

% overlay values
for i = 1:nCond
    for j = 1:nCond
        text(j,i, sprintf('%.2f', C(i,j)), ...
            'Color','w','FontSize',12, ...
            'HorizontalAlignment','center');
    end
end

end
