%% Select recordings that match parameters based on metadata

strainList = {'TU2769'};
internalList = {'IC2','IC6'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 250]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = [1 5]; % frequency of low-pass filter for stimulus command
includeFlag = 1; 

wtCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', includeFlag);

strainList = {'GN381'};

fatCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', includeFlag);


clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep excludeCells;


%% Visual/manual exclusion of bad sweeps 
% This is mainly based on excluding sweeps with leak > 10pA, but also
% sweeps where the recording was lost partway through or some unexpected
% source of noise was clearly at play.

protList ={'WC_Probe8'};
matchType = 'full';

ExcludeSweeps(ephysData, protList, wtCells, 'matchType', matchType);
ExcludeSweeps(ephysData, protList, fatCells, 'matchType', matchType);


%% Find MRCs 
protList ={'WC_Probe8'};
matchType = 'full';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

wtMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'saveSweeps',1,'sweepStats',1);

fatMRCs = IdAnalysis(ephysData,protList,fatCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'saveSweeps',1, 'sweepStats',1);

clear protList sortSweeps matchType

%% Pull out currents vs. stim repetition number

whichMRCs = wtMRCs;
peakCol = 6;
allOns = cell(size(whichMRCs,1),1);
allOffs = cell(size(whichMRCs,1),1);

for iCell = 1:size(whichMRCs,1)
   cellName = whichMRCs{iCell,1};
   allOns{iCell} = whichMRCs{iCell, 4}(:,peakCol);
   allOffs{iCell} = whichMRCs{iCell, 5}(:,peakCol);
end

repLengths = cellfun('length',allOns);
maxLength = max(repLengths);
allOns = cellfun(@(x)vertcat(x,NaN(maxLength-length(x),1)),allOns,'UniformOutput',false);
allOns = cell2mat(allOns');
allOns(allOns==0) = nan;

allOffs = cellfun(@(x)vertcat(x,NaN(maxLength-length(x),1)),allOffs,'UniformOutput',false);
allOffs = cell2mat(allOffs');
allOffs(allOffs==0) = nan;

wtOnNorm = bsxfun(@rdivide,allOns,allOns(1,:));
wtOffNorm = bsxfun(@rdivide,allOffs,allOffs(1,:));

wtRatios = allOffs./allOns;

wtOnMean = nanmean(wtOnNorm,2);
wtOffMean = nanmean(wtOffNorm,2);
wtReps = sum(~isnan(allOns),2);

% Now for fat worms

whichMRCs = fatMRCs;
allOns = cell(size(whichMRCs,1),1);
allOffs = cell(size(whichMRCs,1),1);

for iCell = 1:size(whichMRCs,1)
   cellName = whichMRCs{iCell,1};
   allOns{iCell} = whichMRCs{iCell, 4}(:,peakCol);
   allOffs{iCell} = whichMRCs{iCell, 5}(:,peakCol);
end

repLengths = cellfun('length',allOns);
maxLength = max(repLengths);
allOns = cellfun(@(x)vertcat(x,NaN(maxLength-length(x),1)),allOns,'UniformOutput',false);
allOns = cell2mat(allOns');
allOns(allOns==0) = nan;

allOffs = cellfun(@(x)vertcat(x,NaN(maxLength-length(x),1)),allOffs,'UniformOutput',false);
allOffs = cell2mat(allOffs');
allOffs(allOffs==0) = nan;

fatOnNorm = bsxfun(@rdivide,allOns,allOns(1,:));
fatOffNorm = bsxfun(@rdivide,allOffs,allOffs(1,:));

fatRatios = allOffs./allOns;

fatOnMean = nanmean(fatOnNorm,2);
fatOffMean = nanmean(fatOffNorm,2);
fatReps = sum(~isnan(allOns),2);

%% Plot them

fh(1) = figure(); hold on; ylim(gca,[0 1.75]);
plot(wtOnMean(1:40),'b','LineWidth',2);
plot(fatOnMean(1:40),'Color',[0 1 0],'LineWidth',2);
axh(1) = gca;
fh(2) = figure(); hold on; ylim(gca,[0 1.75]);
plot(wtOffMean(1:40),'b','LineWidth',2);
plot(fatOffMean(1:40),'Color',[0 1 0],'LineWidth',2);
axh(2) = gca;
plotfixer;

plot(axh(1),wtOnNorm(1:40,:),'Color',[.7 .7 1])
plot(axh(1),fatOnNorm(1:40,:),'Color',[.2 .7 .2])
ylabel(axh(1),'Normalized on current');
xlabel(axh(1),'Stimulus #');
chH = get(axh(1),'children');
set(axh(1),'children',flipud(chH));

plot(axh(2),wtOnNorm(1:40,:),'Color',[.7 .7 1])
plot(axh(2),fatOnNorm(1:40,:),'Color',[.2 .7 .2])
ylabel(axh(2),'Normalized off current');
xlabel(axh(2),'Stimulus #');
chH = get(axh(2),'children');
set(axh(2),'children',flipud(chH));


