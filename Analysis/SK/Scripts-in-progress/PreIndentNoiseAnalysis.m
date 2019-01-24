%PreIndentNoiseAnalysis.m


strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'posterior'};

wormPrep = {'dissected'};
cellDist = [40 100];
resistCutoff = '<250';
extFilterFreq = [2.5 5];

noisePreCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq);

% noisePreCells = {'FAT231','FAT232','FAT233','FAT234','FAT235'};

%% Exclusion 
protList ={'NoisePre','WC_Probe8','WC_Probe4','WC_Probe3','WC_Probe12'};
% protList ={'WC_Probe8','WC_Probe4','WC_Probe3','WC_Probe12'};

matchType = 'first';

ExcludeSweeps(ephysData, protList, noisePreCells,'matchType',matchType);
%TODO: Upgrade ExcludeSweeps to allow passing in an existing list, then
%check recording name/series number against it and skip (if overWriteFlag)
%if it's already been vetted. Add new recordings to the end of the list,
%sortrows (and maybe check unique/warn which rows have multiple) and save
%as new file.

%% Running analysis

noisePre = NonStatNoiseAnalysis(ephysData,protList,noisePreCells,...
    'matchType',matchType,'recParameters',ephysMetaData);
clear protList matchType

%% Turn data into tables for Igor

%Get protocol name parts and means/variances for each recording
preRecs = fieldnames(noisePre);
recNames = cell(0); % name of cell
protNames = cell(0); % name of protocol
stimSize = cell(0); % relative displacement of stim
stimPos = cell(0); % total displacement of stim
totMeans = cell(0);
totVars = cell(0);

for iRec = 1:length(preRecs);
    %get rid of empty rows, which mess with sorting strings later
    clearEmpty = arrayfun(@(s) isempty(s.protocol),noisePre.(preRecs{iRec}));
    noisePre.(preRecs{iRec})(clearEmpty)=[];
    
    thesePos = num2cell([noisePre.(preRecs{iRec})(:).position]);
    theseSize = num2cell([noisePre.(preRecs{iRec})(:).size]);

    stimPos = [stimPos thesePos];
    stimSize = [stimSize theseSize];
    protNames = [protNames vertcat({noisePre.(preRecs{iRec})(:).protocol})];
    totMeans = [totMeans vertcat({noisePre.(preRecs{iRec})(:).totalMean})];
    totVars = [totVars vertcat({noisePre.(preRecs{iRec})(:).totalVar})];
    recNames = [recNames repmat(preRecs(iRec),1,length(noisePre.(preRecs{iRec})))];
end

protNames = cellfun(@(x) regexprep(x,'WC_Probe\d*','step'),protNames,'un',0);
protNames = cellfun(@(x,y) regexprep(x,'Noise_PrePulse','pre'),protNames,'un',0);
stimSign = cellfun(@(x) sign(x), stimSize);
stimSign(stimSign<0)=0;
stimSign = num2cell(stimSign);

stimCombine = cellfun(@(x,y,z) sprintf('%dum',x*z + y*~z*-1), stimPos, stimSize, stimSign,'un',0);

stimSign = cellfun(@(x) regexprep(num2str(x),'1','on'),stimSign,'un',0);
stimSign = cellfun(@(x) regexprep(num2str(x),'0','off'),stimSign,'un',0);

[protNames, protIdx] = sort(protNames);
totMeans = totMeans(protIdx);
totVars = totVars(protIdx);
recNames = recNames(protIdx);
stimPos = stimPos(protIdx);
stimCombine = stimCombine(protIdx);
stimSign = stimSign(protIdx);

% Pad with nans and turn means/vars into arrays
meanLengths = cellfun('length',totMeans);
maxLength = max(meanLengths);
totMeans = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totMeans,'UniformOutput',false);
totMeans = cell2mat(totMeans);totVars = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totVars,'UniformOutput',false);
totVars = cellfun(@(x)cat(1,x,NaN(maxLength-length(x),1)),totVars,'UniformOutput',false);
totVars = cell2mat(totVars);

% Create identifying wave names for Igor
waveMeanNames = cellfun(@(w,x,y,z) sprintf('mean_%s_%s_%s_%s', w,x, y,z), protNames, stimCombine, stimSign, recNames, 'un',0);
waveVarNames = cellfun(@(w,x,y,z) sprintf('var_%s_%s_%s_%s', w,x, y,z), protNames, stimCombine, stimSign, recNames, 'un',0);


% write into Excel file for loading into Igor
xlswrite('PatchData/testPre.xls',waveMeanNames,'means');
xlswrite('PatchData/testPre.xls',totMeans,'means','A2');

xlswrite('PatchData/testPre.xls',waveVarNames,'vars');
xlswrite('PatchData/testPre.xls',totVars,'vars','A2');