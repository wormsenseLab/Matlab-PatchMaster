%% Select recordings that match parameters based on metadata

strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = 2.5; % frequency of low-pass filter for stimulus command
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

% This results in a file with three columns: cell ID, series, and sweeps 
% that have been approved for analysis use.

%% Exclusion for current clamp
protList ={'_CC'};
matchType = 'last';

ExcludeSweeps(ephysData, protList, wtCells, 'matchType', matchType);
ExcludeSweeps(ephysData, protList, fatCells, 'matchType', matchType);


%% Find MRPs in current clamp
% antCC_allSteps_wt_2p5(181129).xlsx and fat from 181129.
protList ={'_cc'};
matchType = 'last';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

wtMRPs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'saveZero',1);

fatMRPs = IdAnalysis(ephysData,protList,fatCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'saveZero',1);

clear protList sortSweeps matchType


%% Export for Igor fitting of Boltzmann to each recording

% Set the filename and parameters
fname = 'PatchData/steps_CC_wtVsFat_2p5kHz_allDist(181129).xlsx';

eachSize = [0.5 1 1.5 3 4 5 6 7 8 9 10 11 12]';
distLimits = [40 200]; % limit to same average distance for anterior and posterior
distCol = 12;
whichStim = 1; %on or off

% Step to normalize to
normFlag = 1; % yes, add normalized sheets
normStepSize = 10; % step size to normalize to
corrOn = 0; %no voltage attenuation correction if 0

whichMRPs = wtMRPs;
dType = 'decay'; 

switch dType
    case 'pot'
        peakCol = 6; % peak potential
    case 'act'
        peakCol = 8; % activation time constant
    case 'decay'
        peakCol = 9; % decay time constant
end

if whichStim == 2
    eachSize = -eachSize;
    normStepSize = -normStepSize;
end

thisAtt = attenuationData(:,[2 8 10]);
thisName = cell(0);
wt_Out = cell(length(eachSize)+1,size(whichMRPs,1));

for iCell = 1:size(whichMRPs,1)
    thisCell = whichMRPs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRPs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        wt_Out{1,iCell} = [dType '_' thisName{iCell,1}];
        hasAtt = strcmp(thisName{iCell,1},thisAtt(:,1));
        Iact = [];
        
        for iSize = 1:length(eachSize)
            stepSize = eachSize(iSize);
            whichStep = round(thisCell(:,1)*2)/2 == stepSize; %round to nearest 0.5
            if any(whichStep)
                
                Im = thisCell(whichStep,peakCol);
                Iact(iSize,1) = Im;
                
            else
                Iact(iSize,1) = nan;
            end
        end
        wt_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end
wt_Name = [{'wt', 'wt_Dist'}; thisName];
wt_Out = wt_Out(:,~cellfun(@isempty, wt_Out(1,:))); % clear out empty waves (where dist didn't match)
wt_Name = wt_Name(~cellfun(@isempty, wt_Name(:,1)),:);

wt_Out = [[{'stepSize'};num2cell(eachSize)] wt_Out]; % append stepSize wave

thisName = cell(0);
whichMRPs = fatMRPs;
fat_Out = cell(length(eachSize)+1,size(whichMRPs,1));

for iCell = 1:size(whichMRPs,1)
    thisCell = whichMRPs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRPs{iCell,1}; %name
        thisName{iCell,2} = thisDist; %distance
        fat_Out{1,iCell} = [dType '_' thisName{iCell,1}];
        Iact = [];
        
        for iSize = 1:length(eachSize)
            stepSize = eachSize(iSize);
            whichStep = round(thisCell(:,1)*2)/2 == stepSize; %round to nearest 0.5
            if any(whichStep)
                Iact(iSize,1) = thisCell(whichStep,peakCol);
            else
                Iact(iSize,1) = nan;
            end
        end
        fat_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end

fat_Name = [{'fat', 'fat_Dist'}; thisName];
fat_Out = fat_Out(:,~cellfun(@isempty, fat_Out(1,:)));
fat_Name = fat_Name(~cellfun(@isempty, fat_Name(:,1)),:);

fat_Out = [[{'stepSize'};num2cell(eachSize)] fat_Out]; % append stepSize wave


xlswrite(fname,wt_Out,['wt_' dType]);
xlswrite(fname,fat_Out,['fat_' dType]);
xlswrite(fname,wt_Name,'wtStats');
xlswrite(fname,fat_Name,'fatStats');


clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName


% Normalize based on current/charge at given step size and write to xls
if normFlag && strcmp(dType,'pot')
    %WT
    which_Out = wt_Out;
    norm_Out = cell(size(which_Out));
    
    currMat = cell2mat(which_Out(2:end,2:end));
    whichRow = find(cellfun(@(x) x==normStepSize, which_Out(2:end,1)));
    normRow = currMat(whichRow,:);
    normMat = bsxfun(@rdivide,currMat,normRow);
    norm_Out(2:end,2:end) = num2cell(normMat);
    norm_Out(:,1) = which_Out(:,1);
    norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);
    
    
    wt_Norm = norm_Out;
    wt_Max = normRow';
    
    % FAT
    which_Out = fat_Out;
    norm_Out = cell(size(which_Out));
    
    currMat = cell2mat(which_Out(2:end,2:end));
    whichRow = find(cellfun(@(x) x==normStepSize, which_Out(2:end,1)));
    normRow = currMat(whichRow,:);
    normMat = bsxfun(@rdivide,currMat,normRow);
    norm_Out(2:end,2:end) = num2cell(normMat);
    norm_Out(:,1) = which_Out(:,1);
    norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);
    
    fat_Norm = norm_Out;
    fat_Max = normRow';
    
    xlswrite(fname,wt_Norm,['wt_' dType '_Norm']);
    xlswrite(fname,fat_Norm,['fat_' dType '_Norm']);
    
    clear currMat which Row normRow normMat norm_Out
end

%% Pull out and plot sweeps at one displacement

whichMRPs = wtMRPs;
stepSize = 5;

for iCell = 1:size(whichMRPs,1)
    thisCell = whichMRPs{iCell,whichStim+2};
    whichStep = round(thisCell(:,1)*2)/2 == stepSize; %round to nearest 0.5
    thisTrace = whichMRPs{iCell,2}(whichStep,:);
end


%% Pull out capacitances for comparison

whichMRPs = wtMRPs;
thisCap = [];
for iCell = 1:size(whichMRPs,1)    
    thisName = whichMRPs{iCell,1}; %name
    recRow = strcmp(thisName,ephysMetaData(:,1));
    capCol = strcmp('C_in (pF)',ephysMetaData(1,:));
    thisCap(iCell,1) = ephysMetaData{recRow,capCol};
end

wtCap = thisCap;

whichMRPs = fatMRPs;
thisCap = [];
for iCell = 1:size(whichMRPs,1)    
    thisName = whichMRPs{iCell,1}; %name
    recRow = strcmp(thisName,ephysMetaData(:,1));
    capCol = strcmp('C_in (pF)',ephysMetaData(1,:));
    thisCap(iCell,1) = ephysMetaData{recRow,capCol};
end

fatCap = thisCap;
