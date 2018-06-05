% 8um Noise Analysis

protList ={'WC_Probe1','WC_Probe4','wC_Probe3','WC_Probe8'};
matchType = 'full';
strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'posterior'};
wormPrep = {'dissected'};
cellDist = [40 250];
resistCutoff = '<250';
extFilterFreq = 2.5;

posteriorNoiseCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
     'stimFilterFrequencykHz', extFilterFreq);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%%

%NEXT: run exclude sweeps and select sweeps
%THEN: run NonStatNoiseAnalysis.m
%THEN: take 8um 
%LATER: integrate thisDist into Id and NonStat as something you can sort by

ExcludeSweeps(ephysData, protList, posteriorNoiseCells, 'matchType', matchType);

clear protList matchType;

%%
protList ={'WC_Probe1','WC_Probe4','wC_Probe3','WC_Probe8'};
matchType = 'full';
noisePost = NonStatNoiseAnalysis(ephysData,protList,posteriorNoiseCells,'matchType',matchType,...
    'recParameters',ephysMetaData);
clear protList matchType

%% Turn data into tables for Igor

%Get protocol name parts and means/variances for each recording
noiseRecs = fieldnames(noisePost);
recNames = cell(0); % name of cell
protNames = cell(0); % name of protocol
stimNames = cell(0); % number of stim (e.g., 1=on, 2=off) as string
totMeans = cell(0);
totVars = cell(0);

for iRec = 1:length(noiseRecs);
    %get rid of empty rows, which mess with sorting strings later
    clearEmpty = arrayfun(@(s) isempty(s.protocol),noisePost.(noiseRecs{iRec}));
    noisePost.(noiseRecs{iRec})(clearEmpty)=[];
    
    theseStim = cellfun(@(x) sprintf('stim%d',x),...
        num2cell([noisePost.(noiseRecs{iRec})(:).stimNum]),'un',0);

    stimNames = [stimNames theseStim];
    protNames = [protNames vertcat({noisePost.(noiseRecs{iRec})(:).protocol})];
    totMeans = [totMeans vertcat({noisePost.(noiseRecs{iRec})(:).totalMean})];
    totVars = [totVars vertcat({noisePost.(noiseRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(noiseRecs(iRec),1,length(noisePost.(noiseRecs{iRec})))];
end

[protNames, protIdx] = sort(protNames);
totMeans = totMeans(protIdx);
totVars = totVars(protIdx);
recNames = recNames(protIdx);
stimNames = stimNames(protIdx);

% Pad with nans and turn means/vars into arrays
meanLengths = cellfun('length',totMeans);
maxLength = max(meanLengths);
totMeans = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totMeans,'UniformOutput',false);
totMeans = cell2mat(totMeans);
totVars = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totVars,'UniformOutput',false);
totVars = cell2mat(totVars);

protTimes = cellfun(@(x) sprintf('%sms',x(isstrprop(x,'digit'))), protNames,'un',0);

% Create identifying wave names for Igor
waveMeanNames = cellfun(@(x,y,z) sprintf('mean_%s_%s_%s', x, y,z), protTimes, stimNames, recNames, 'un',0);
waveVarNames = cellfun(@(x,y,z) sprintf('var_%s_%s_%s', x, y,z), protTimes, stimNames, recNames, 'un',0);


% write into Excel file for loading into Igor
xlswrite('PatchData/noisePost.xls',waveMeanNames,'means');
xlswrite('PatchData/noisePost.xls',totMeans,'means','A2');

xlswrite('PatchData/noisePost.xls',waveVarNames,'vars');
xlswrite('PatchData/noisePost.xls',totVars,'vars','A2');

