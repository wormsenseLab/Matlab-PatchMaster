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

%Get protocol names and means/variances for each recording
trapRecs = fieldnames(noiseTrap);
recNames = cell(0);
protNames = cell(0);
totMeans = cell(0);
totVars = cell(0);
for iRec = 1:length(trapRecs);
    protNames = [protNames vertcat({noiseTrap.(trapRecs{iRec})(:).protocol})];
    totMeans = [totMeans vertcat({noiseTrap.(trapRecs{iRec})(:).totalMean})];
    totVars = [totMeans vertcat({noiseTrap.(trapRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(trapRecs(iRec),1,length(noiseTrap.(trapRecs{iRec})))];
end

% Deal with empty cells so you can sort as an array of strings
tf = cellfun('isempty',protNames);
protNames(tf) = {''};

[protNames, protIdx] = sort(protNames);
totMeans = totMeans(protIdx);
totVars = totVars(protIdx);
recNames = recNames(protIdx);

% Pad with nans and turn means/vars into arrays
meanLengths = cellfun('length',totMeans);
maxLength = max(meanLengths);
totMeans = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totMeans,'UniformOutput',false);
totMeans = cell2mat(totMeans);
totVars = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totVars,'UniformOutput',false);
totVars = cell2mat(totVars);

% Drop the empty columns
totMeans = totMeans(~cellfun('isempty',protNames));
totVars = totVars(~cellfun('isempty',protNames));
recNames = recNames(~cellfun('isempty',protNames));
protNames = protNames(~cellfun('isempty',protNames));

protTimes = cellfun(@(x) sprintf('%sms',x(isstrprop(x,'digit'))), protNames,'un',0);

% Create identifying wave names for Igor
waveMeanNames = cellfun(@(x,y) sprintf('mean_%s_%s', x, y), protTimes, recNames, 'un',0);
waveVarNames = cellfun(@(x,y) sprintf('var_%s_%s', x, y), protTimes, recNames, 'un',0);


% TODO: figure out how to add header names for the waves
fid = fopen('PatchData/testTrap.txt','w',
dlmwrite('PatchData/testTrap.txt',waveMeanNames,'delimiter','\t'); 