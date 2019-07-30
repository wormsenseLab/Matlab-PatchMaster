% SpeedControlledSteps.m

%%Import data and metadata
ephysData = ImportPatchData('incl',1);
projects = {'FAT';'SYM'};

%fix a couple messed up names
ephysData.FAT029 = ephysData.FAT029s;
ephysData.FAT017 = ephysData.FAT017e001;
ephysData.FAT164 = ephysData.FAT164001;
ephysData = rmfield(ephysData,{'FAT017e001';'FAT029s';'FAT164001'});

ephysData = FilterProjectData(ephysData, projects);
clear projects;

ephysMetaData = ImportMetaData();  %Recording Database.xlsx
attenuationData = ImportMetaData(); %AttenuationCalcs.xlsx

%% Select recordings that match parameters based on metadata

strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = 2.5; % frequency of low-pass filter for stimulus command
includeFlag = 1; 

anteriorControlCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', includeFlag);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep excludeCells;

% This results in a file with three columns: cell ID, series, and sweeps 
% that have been approved for analysis use.

%% Visual/manual exclusion of bad sweeps 
% This is mainly based on excluding sweeps with leak > 10pA, but also
% sweeps where the recording was lost partway through or some unexpected
% source of noise was clearly at play.

% If you have already made a list of sweeps, skip this section and load the
% list in the next section.

protList ={'Control'};
matchType = 'first';

ExcludeSweeps(ephysData, protList, anteriorControlCells, 'matchType', matchType);


%% Find MRCs 
% Lists: antAllSteps.xlsx and postAllSteps.xlsx from 180810 should contain
% all relevant steps for anterior and posterior stimulation, respectively.
protList ={'Control'};
matchType = 'first';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

[anteriorMRCs, antStim] = IdAnalysis(ephysData,protList,anteriorControlCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

clear protList sortSweeps matchType

%% Average and plot

% Set the filename
fname = 'PatchData/controlSpeed_attCorrected(190729).xls';

eachSize = [4 8 12]';
distLimits = [40 200]; % limit to same average distance for anterior and posterior
distCol = 12;
whichStim = 1; %on or off
noCorr = 0; % 1 to skip space clamp voltage attenuation correction
Vc = -0.06; %in V
Ena = 0.094; % in V

% Manually change dType to 'curr', 'char', 'act', or 'decay' depending on
% which analysis you would like to do.
whichMRCs = anteriorMRCs;
dType = 'curr'; 

switch dType
    case 'curr'
        peakCol = 6; % peak current
    case 'char'
        peakCol = 11; % integrated current/charge
    case 'act'
        peakCol = 8; % activation time constant
    case 'decay'
        peakCol = 9; % decay time constant
end

if whichStim == 2
    eachSize = -eachSize;
end

thisAtt = attenuationData(:,[2 8 10]);
thisName = cell(0);
ant_Out = cell(length(eachSize)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        ant_Out{1,iCell} = [dType '_' thisName{iCell,1}];
        hasAtt = strcmp(thisName{iCell,1},thisAtt(:,1));
        Iact = [];
        
        for iSize = 1:length(eachSize)
            stepSize = eachSize(iSize);
            whichStep = round(thisCell(:,1)*2)/2 == stepSize; %round to nearest 0.5
            if any(whichStep)
                
                if any(hasAtt) && thisAtt{hasAtt,2} %if attenuation calc exists and not omitCell
                    if noCorr == 0 || strcmp(dType,'peak') || strcmp(dType,'charge') %attenuation correction for current/charge but not taus
                        Vm = Vc * thisAtt{hasAtt,3};
                        Im = thisCell(whichStep,peakCol);
                       
                        Iact(iSize,1) = (Im * (Vc-Ena)) ./ (Vm-Ena);
                    else
                        Im = thisCell(whichStep,peakCol);
                        Iact(iSize,1) = Im;
                    end
                else
                    Iact(iSize,1) = nan;
                end
       
            else
                Iact(iSize,1) = nan;
            end
        end
        ant_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end
ant_Name = [{'Anterior', 'Ant_Dist'}; thisName];
ant_Out = ant_Out(:,~cellfun(@isempty, ant_Out(1,:))); % clear out empty waves (where dist didn't match)
ant_Name = ant_Name(~cellfun(@isempty, ant_Name(:,1)),:);

ant_Out = [[{'stepSize'};num2cell(eachSize)] ant_Out]; % append stepSize wave

xlswrite(fname,ant_Out,['ant_' dType]);
xlswrite(fname,ant_Name,'antStats');

clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName

%% Normalize based on current/charge at given step size and write to xls
% Normalizes within each condition

stepSize = 12; % must be negative for off

%Anterior
which_Out = ant_Out;
norm_Out = cell(size(which_Out));

currMat = cell2mat(which_Out(2:end,2:end));
whichRow = find(cellfun(@(x) abs(x - stepSize) < 0.5, which_Out(2:end,1)));
normRow = currMat(whichRow,:);
normMat = bsxfun(@rdivide,currMat,normRow);
norm_Out(2:end,2:end) = num2cell(normMat);
norm_Out(:,1) = which_Out(:,1);
norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);

ant_Norm = norm_Out;
xlswrite(fname,ant_Norm,['ant_' dType '_Norm']);

clear currMat which Row normRow normMat norm_Out


