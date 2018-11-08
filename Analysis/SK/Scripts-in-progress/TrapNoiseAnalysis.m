%TrapNoiseAnalysis.m


strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 100];
resistCutoff = '<250';
extFilterFreq = [2.5 5];

noiseTrapCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', 1);

%% Exclusion 
protList ={'WCProbe1','WCProbe3'};

matchType = 'first';

ExcludeSweeps(ephysData, protList, noiseTrapCells,'matchType',matchType);
%TODO: Upgrade ExcludeSweeps to allow passing in an existing list, then
%check recording name/series number against it and skip (if overWriteFlag)
%if it's already been vetted. Add new recordings to the end of the list,
%sortrows (and maybe check unique/warn which rows have multiple) and save
%as new file.

%% Approving Traces (old version)
% protList ={'Noise_Trap','WCProbe8'};
% matchType = 'first';
% strainList = {'TU2769'};
% internalList = {'IC6'};
% stimPosition = {'posterior'};
% wormPrep = {'dissected'};
% 
% % define newCells when importing additional data so you can just add to the
% % list of approved sweeps
% 
% trapCells = FilterRecordings(ephysData, ephysMetaData, newCells, ...
%     'strain', strainList, 'internal', internalList, ...
%      'stimLocation', stimPosition, 'wormPrep', wormPrep);
% 
% ExcludeSweeps(ephysData, protList, trapCells, 'matchType', matchType);
% 
% clear protList strainList internalList cellTypeList stimPosition matchType ans wormPrep;

%% Running analysis

protList ={'NoiseTrap','WC_Probe8'};
matchType = 'first';
noiseTrap = NonStatNoiseAnalysis(ephysData,protList,noiseTrapCells,'matchType',matchType,...
    'recParameters',ephysMetaData);
clear protList matchType

%% Turn data into tables for Igor

fname = 'PatchData/noiseTrapAnt(181001).xls';

%Get protocol name parts and means/variances for each recording
trapRecs = fieldnames(noiseTrap);
recNames = cell(0); % name of cell
protVel = []; % name of protocol
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
    protVel = [protVel; vertcat(noiseTrap.(trapRecs{iRec})(:).velocity)];
    totMeans = [totMeans vertcat({noiseTrap.(trapRecs{iRec})(:).totalMean})];
    totVars = [totVars vertcat({noiseTrap.(trapRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(trapRecs(iRec),1,length(noiseTrap.(trapRecs{iRec})))];
end

[protVel, protIdx] = sort(protVel);
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
protVel = num2cell(abs(protVel))';

% Create identifying wave names for Igor
waveMeanNames = cellfun(@(x,y,z) sprintf('mean_%s_%s_%dums', x, y,z), recNames, stimNames, protVel, 'un',0);
waveVarNames = cellfun(@(x,y,z) sprintf('var_%s_%s_%dums', x, y,z), recNames, stimNames, protVel, 'un',0);


% write into Excel file for loading into Igor
xlswrite(fname,waveMeanNames,'means');
xlswrite(fname,totMeans,'means','A2');

xlswrite(fname,waveVarNames,'vars');
xlswrite(fname,totVars,'vars','A2');