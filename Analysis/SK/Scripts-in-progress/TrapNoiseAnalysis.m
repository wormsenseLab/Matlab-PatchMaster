%TrapNoiseAnalysis.m

%% Approving Traces
protList ={'Noise_Trap'};
matchType = 'first';
strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'posterior'};
wormPrep = {'dissected'};

% define newCells when importing additional data so you can just add to the
% list of approved sweeps

trapCells = FilterRecordings(ephysData, ephysMetaData, newCells, ...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep);

ExcludeSweeps(ephysData, protList, trapCells, 'matchType', matchType);

clear protList strainList internalList cellTypeList stimPosition matchType ans wormPrep;

%% Running analysis

protList ={'NoiseTrap'};
matchType = 'first';
noiseTrap = NonStatNoiseAnalysis(ephysData,protList,trapCells,'matchType',matchType);
clear protList matchType

%% Turn data into tables for Igor

%Get protocol name parts and means/variances for each recording
trapRecs = fieldnames(noiseTrap);
recNames = cell(0); % name of cell
protNames = cell(0); % name of protocol
stimNames = cell(0); % number of stim (e.g., 1=on, 2=off) as string
totMeans = cell(0);
totVars = cell(0);

for iRec = 1:length(trapRecs);
    %get rid of empty rows, which mess with sorting strings later
    clearEmpty = arrayfun(@(s) isempty(s.protocol),noiseTrap.(trapRecs{iRec}));
    noiseTrap.(trapRecs{iRec})(clearEmpty)=[];
    
    theseStim = cellfun(@(x) sprintf('stim%d',x),...
        num2cell([noiseTrap.(trapRecs{iRec})(:).stimNum]),'un',0);

    stimNames = [stimNames theseStim];
    protNames = [protNames vertcat({noiseTrap.(trapRecs{iRec})(:).protocol})];
    totMeans = [totMeans vertcat({noiseTrap.(trapRecs{iRec})(:).totalMean})];
    totVars = [totVars vertcat({noiseTrap.(trapRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(trapRecs(iRec),1,length(noiseTrap.(trapRecs{iRec})))];
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
xlswrite('PatchData/testTrap.xls',waveMeanNames,'means');
xlswrite('PatchData/testTrap.xls',totMeans,'means','A2');

xlswrite('PatchData/testTrap.xls',waveVarNames,'vars');
xlswrite('PatchData/testTrap.xls',totVars,'vars','A2');