% ChargeExploration.m
% 
% 
% Exploring asymmetry in charge with multi-rep noise analysis traces.
% Variation in charge asymmetry with small steps, and whether charge
% actually returns to zero by the end of the on step even if there is no
% current detectable above noise.

strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200];
resistCutoff = '<250';
extFilterFreq = [2.5 5];

smallStepCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', 1);

clear strainList internalList stimPosition wormPrep cellDist resistCutoff extFilterFreq
%% Exclusion 
protList ={'WCProbe1','WCProbe3','WCProbe8'};

matchType = 'full';

ExcludeSweeps(ephysData, protList, smallStepCells,'matchType',matchType);

%% Running noise analysis

protList ={'WCProbe1'};
matchType = 'full';
steps1 = NonStatNoiseAnalysis(ephysData,protList,smallStepCells,'matchType',matchType,...
    'recParameters',ephysMetaData,'responseTime',250);

protList ={'WCProbe3'};
steps3 = NonStatNoiseAnalysis(ephysData,protList,smallStepCells,'matchType',matchType,...
    'recParameters',ephysMetaData,'responseTime',250);


protList ={'WCProbe8'};
steps8 = NonStatNoiseAnalysis(ephysData,protList,smallStepCells,'matchType',matchType,...
    'recParameters',ephysMetaData,'responseTime',250);


clear protList matchType


%% Calculate on/off charge for each trace of each recording
whichSteps = steps1;
theseNames = fieldnames(whichSteps);
peaks = cell(0);
charges = cell(0);
cellPeaks = cell(0);
peakRatios = cell(0);
chargeRatios = cell(0);
sf = 5; %khz
dataType = 'A';
tauType = 'thalfmax';
integrateFlag = 1;

for iCell = 1:length(theseNames)
    cellName = theseNames{iCell};
    for iSeries = 1:2
        theseTraces(:,:,iSeries) = whichSteps.(cellName)(iSeries).traces';
        stimParams = [whichSteps.(cellName)(iSeries).stimLoc, whichSteps.(cellName)(iSeries).size, ...
            whichSteps.(cellName)(iSeries).size, whichSteps.(cellName)(iSeries).position, ...
            whichSteps.(cellName)(iSeries).velocity, 1, whichSteps.(cellName)(iSeries).distance];
        stimParams = repmat(stimParams,size(theseTraces,1),1);
        
        cellPeaks{iCell,iSeries} = findMRCs(stimParams, theseTraces(:,:,iSeries), sf, dataType, ...
            'tauType', tauType, 'integrateCurrent',integrateFlag,'threshTime',15);
    end
    peaks{iCell}(:,1) = cellPeaks{iCell,1}(:,6);
    peaks{iCell}(:,2) = cellPeaks{iCell,2}(:,6);
    charges{iCell}(:,1) = cellPeaks{iCell,1}(:,7);
    charges{iCell}(:,2) = cellPeaks{iCell,2}(:,7);
    peakRatios{iCell} = cellPeaks{iCell,2}(:,6)./cellPeaks{iCell,1}(:,6);
    chargeRatios{iCell} = cellPeaks{iCell,2}(:,7)./cellPeaks{iCell,1}(:,7);
end

% Pad and concatenate for Igor assuming max length is 64 sweeps
for i = 1:length(peakRatios)
   peakRatios{i} = vertcat(peakRatios{i}, nan(64-length(peakRatios{i}),1));
   chargeRatios{i} = vertcat(chargeRatios{i}, nan(64-length(chargeRatios{i}),1));
   peaks{i} = vertcat(peaks{i}, nan(64-length(peaks{i}),2));
   charges{i} = vertcat(charges{i}, nan(64-length(charges{i}),2));
end
peakRatioCat = [peakRatios{:}]';
chargeRatioCat = [chargeRatios{:}]';

peakCat = [];
chargeCat = [];
for i = 1:2
   a = [peaks{:}];
   peakCat(:,:,i) = a(:,i:2:end);
   a = [charges{:}];
   chargeCat(:,:,i) = a(:,i:2:end);
end

%%
figure();
for i = 1:2
    subplot(1,2,i);
    plot((theseTraces(:,:,iSeries)+repmat(offset',1,1328))')
    vline(76)
    xlim([0 1000]);
    ylim([-70e-11,.5e-11]); 
end