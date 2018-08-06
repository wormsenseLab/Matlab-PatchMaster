%PreIndentNoiseAnalysis.m


strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 100];
resistCutoff = '<251';
extFilterFreq = [2.5 5];

noisePreCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq);

noisePreCells = {'FAT231','FAT232','FAT233','FAT234','FAT235'};

%% Exclusion 
protList ={'NoisePre','WC_Probe'};

matchType = 'first';

ExcludeSweeps(ephysData, protList, testCell,'matchType',matchType);
%TODO: Upgrade ExcludeSweeps to allow passing in an existing list, then
%check recording name/series number against it and skip (if overWriteFlag)
%if it's already been vetted. Add new recordings to the end of the list,
%sortrows (and maybe check unique/warn which rows have multiple) and save
%as new file.

%% Running analysis

noisePre = NonStatNoiseAnalysis(ephysData,protList,noisePreCells,'matchType',matchType);
clear protList matchType

%% Turn data into tables for Igor

%Get protocol name parts and means/variances for each recording
preRecs = fieldnames(noisePre);
recNames = cell(0); % name of cell
protNames = cell(0); % name of protocol
stimNames = cell(0); % number of stim (e.g., 1=on, 2=off) as string
totMeans = cell(0);
totVars = cell(0);

for iRec = 1:length(preRecs);
    %get rid of empty rows, which mess with sorting strings later
    clearEmpty = arrayfun(@(s) isempty(s.protocol),noisePre.(preRecs{iRec}));
    noisePre.(preRecs{iRec})(clearEmpty)=[];
    
    theseStim = cellfun(@(x) sprintf('stim%d',x),...
        num2cell([noisePre.(preRecs{iRec})(:).stimNum]),'un',0);

    stimNames = [stimNames theseStim];
    protNames = [protNames vertcat({noisePre.(preRecs{iRec})(:).protocol})];
    totMeans = [totMeans vertcat({noisePre.(preRecs{iRec})(:).totalMean})];
    totVars = [totVars vertcat({noisePre.(preRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(preRecs(iRec),1,length(noisePre.(preRecs{iRec})))];
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