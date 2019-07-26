% IntactAnalysis.m 
% I-d and Q-d curves for anterior dissected vs. anterior intact
% preparations

%%Import data and metadata
ephysData = ImportPatchData('incl',1);
projects = {'FAT';'SYM';'SPC'};

%fix a couple messed up names
ephysData.FAT029 = ephysData.FAT029s;
ephysData.FAT017 = ephysData.FAT017e001;
ephysData.FAT164 = ephysData.FAT164001;
ephysData = rmfield(ephysData,{'FAT017e001';'FAT029s';'FAT164001'});

ephysData = FilterProjectData(ephysData, projects);
clear projects;

% Get metadata from sheet
ephysMetaData = ImportMetaData();  %Recording Database.xlsx
attenuationData = ImportMetaData(); %AttenuationCalcs.xlsx

%% Analyze capacity transient for C, Rs, and tau

ephysData = CtAnalysis(ephysData);

clear lastfit;

%% Select recordings that match parameters based on metadata

strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = 2.5; % frequency of low-pass filter for stimulus command
includeFlag = 1; 
% cultT = 15;

dissectedDistCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', includeFlag);



wormPrep = {'intact'};

intactDistCells = FilterRecordings(ephysData, ephysMetaData,...
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

protList ={'WC_Probe';'WC_ProbeSmall';'WC_ProbeLarge';'NoPrePulse'};
matchType = 'full';

ExcludeSweeps(ephysData, protList, dissectedDistCells, 'matchType', matchType);
% ExcludeSweeps(ephysData, protList, intactDistCells, 'matchType', matchType);

%% Find MRCs 
% antAllSteps.xlsx and postAllSteps.xlsx from 180810
protList ={'WC_Probe';'NoPre'};
matchType = 'first';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

[dissectedMRCs, dissectedStim] = IdAnalysis(ephysData,protList,dissectedDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

% Some intact recordings have both WC_Probe and NoPre, which findMRCs
% cannot yet deal with (combineSweeps has not been set up). Different
% timepoints for the same magnitude/speed step cannot be done easily right
% now, so only use WC_Probe steps for intact cells.
protList ={'WC_Probe'}; 

intactMRCs = IdAnalysis(ephysData,protList,intactDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

clear protList sortSweeps matchType



%% Correct all sizes and export for Igor plotting

% Set the filename
fname = 'PatchData/intVsDiss_off_allDist_subQ(190624).xls';


eachSize = [0.5 1 1.5 3 4 5 6 7 8 9 10 11 12]';
distLimits = [40 200]; % limit to same average distance for anterior and posterior
distCol = 12;
whichStim = 2; %on or off
noCorr = 0; % 1 to skip space clamp voltage attenuation correction
Vc = -0.06; %in V
Ena = 0.094; % in V


whichMRCs = dissectedMRCs;
dType = 'decay'; 

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
diss_Out = cell(length(eachSize)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        diss_Out{1,iCell} = [dType '_' thisName{iCell,1}];
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
        diss_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end
diss_Name = [{'Dissected', 'Dissected_Dist'}; thisName];
diss_Out = diss_Out(:,~cellfun(@isempty, diss_Out(1,:))); % clear out empty waves (where dist didn't match)
diss_Name = diss_Name(~cellfun(@isempty, diss_Name(:,1)),:);

diss_Out = [[{'stepSize'};num2cell(eachSize)] diss_Out]; % append stepSize wave

thisName = cell(0);
whichMRCs = intactMRCs;
int_Out = cell(length(eachSize)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        int_Out{1,iCell} = [dType '_' thisName{iCell,1}];
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
        int_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end

int_Name = [{'Intact', 'Intact_Dist'}; thisName];
int_Out = int_Out(:,~cellfun(@isempty, int_Out(1,:)));
int_Name = int_Name(~cellfun(@isempty, int_Name(:,1)),:);

int_Out = [[{'stepSize'};num2cell(eachSize)] int_Out]; % append stepSize wave



xlswrite(fname,diss_Out,['diss_' dType]);
xlswrite(fname,int_Out,['int_' dType]);
xlswrite(fname,diss_Name,'dissStats');
xlswrite(fname,int_Name,'intStats');


clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName


%% Normalize based on current/charge at given step size and write to xls
% Normalizes within each condition

stepSize = -10; % must be negative for off

%Dissected
which_Out = diss_Out;
norm_Out = cell(size(which_Out));

currMat = cell2mat(which_Out(2:end,2:end));
whichRow = find(cellfun(@(x) x==stepSize, which_Out(2:end,1)));
normRow = currMat(whichRow,:);
normMat = bsxfun(@rdivide,currMat,normRow);
norm_Out(2:end,2:end) = num2cell(normMat);
norm_Out(:,1) = which_Out(:,1);
norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);


diss_Norm = norm_Out;


% Intact
which_Out = int_Out;
norm_Out = cell(size(which_Out));

currMat = cell2mat(which_Out(2:end,2:end));
whichRow = find(cellfun(@(x) x==stepSize, which_Out(2:end,1)));
normRow = currMat(whichRow,:);
normMat = bsxfun(@rdivide,currMat,normRow);
norm_Out(2:end,2:end) = num2cell(normMat);
norm_Out(:,1) = which_Out(:,1);
norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);

int_Norm = norm_Out;

xlswrite(fname,diss_Norm,['diss_' dType '_Norm']);
xlswrite(fname,int_Norm,['int_' dType '_Norm']);

clear currMat which Row normRow normMat norm_Out


%%
dissOffChar = diss_Out;
intOffChar = int_Out;

%% Calculate ratios
% re-save ant_Out as antOff/antOnCurr for the corresponding cases
fname = 'PatchData/intVsDiss_ratio_allDist_subQ(190624).xls';
dType = 'curr';

on = intOnCurr;
off = intOffCurr;
out = on(1,:);

onVel = [on{2:end,1}];
offVel = -[off{2:end,1}];


for i = 1:length(onVel)
    out{i+1,1} = onVel(i);
    out(i+1,2:end) = cellfun(@(x,y) x./y, off(i+1,2:end), on(i+1,2:end),'un',0);
end

out = out(~cellfun(@isempty,out(:,1)),:);
out(1,:) = cellfun(@(x) regexprep(x,'stim1','ratio'),out(1,:),'un',0);

xlswrite(fname,out,['int_' dType '_Ratio']);




