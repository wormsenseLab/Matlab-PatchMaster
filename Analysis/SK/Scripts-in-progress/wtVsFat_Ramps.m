% I-d and Q-d curves for anterior (corrected for space clamp) vs. posterior
% Current or charge vs. displacement, and vs. distance 

%%Import data and metadata
ephysData = ImportPatchData('incl',1);
projects = {'FAT';'SYM'};
ephysData = FilterProjectData(ephysData, projects);
clear projects;

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
cultT = 15;
includeFlag = 1; 

wtCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'TCultC', cultT, ...
    'included', includeFlag);

strainList = {'GN381'};

fatCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'TCultC', cultT, ...
    'included', includeFlag);


clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep excludeCells;

% This results in a file with three columns: cell ID, series, and sweeps 
% that have been approved for analysis use.

%% Visual/manual exclusion of bad sweeps 
% This is mainly based on excluding sweeps with leak > 10pA, but also
% sweeps where the recording was lost partway through or some unexpected
% source of noise was clearly at play.

protList ={'DispRate'};
matchType = 'first';

ExcludeSweeps(ephysData, protList, wtCells, 'matchType', matchType);
ExcludeSweeps(ephysData, protList, fatCells, 'matchType', matchType);


% %% Drop super slow ramps that don't even work
% protList ={'DispRate'};
% matchType = 'first';
% 
% ExcludeSweeps(ephysData, protList, wtCells, 'matchType', matchType,'channel',2);
% ExcludeSweeps(ephysData, protList, fatCells, 'matchType', matchType,'channel',2);
% 
% 


%% Find MRCs 
% antAllRampsWT(180810).xlsx and antAllRampsFat(181129).xlsx
protList ={'DispRate'};
matchType = 'first';

sortSweeps = {'velocity','velocity','magnitude','magnitude'};

wtMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

fatMRCs = IdAnalysis(ephysData,protList,fatCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

clear protList sortSweeps matchType



%% Correct all sizes and export for Igor fitting of Boltzmann to each recording

% Set the filename and parameters
fname = 'PatchData/ramps_wtVsFat_allFreq_allDist(190218).xlsx';

tol = 12;

allVels = vertcat(wtMRCs{:,3});
allVels = allVels(:,1);
[~,~,eachVel] = sortRowsTol(allVels,tol,1);
distLimits = [40 200]; % limit to same average distance for anterior and posterior
distCol = 12;
whichStim = 1; %on or off
Vc = -0.06; %in V
Ena = 0.094; % in V
extFiltCorr = 1; % if on, accounts for non-stimCom-filtered recordings with external filter
                 % by correcting everything above 20mm/s to 20mm/s.

% Step to normalize to
normFlag = 1; % yes, add normalized sheets
normVel = 19820; % must be negative for off
corrOn = 0; %no voltage attenuation correction if 0

whichMRCs = wtMRCs;
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
    eachVel = -eachVel;
end
if extFiltCorr
    eachVel(eachVel>19820)=19820;
    eachVel = unique(eachVel);
end

thisAtt = attenuationData(:,[2 8 10]);
thisName = cell(0);
wt_Out = cell(length(eachVel)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        wt_Out{1,iCell} = [dType '_' thisName{iCell,1}];
        hasAtt = strcmp(thisName{iCell,1},thisAtt(:,1));
        if extFiltCorr
           whichCorr=thisCell(:,1)>20000;
           thisCell(whichCorr,1) = 19820;
        end
        Iact = [];
        
        for iSize = 1:length(eachVel)
            stepSize = eachVel(iSize);
            whichStep = abs(roundVel(thisCell(:,1))-stepSize)<tol;
            if any(whichStep)
                
                if any(hasAtt) && thisAtt{hasAtt,2} %if attenuation calc exists and not omitCell
                    if corrOn == 1 && strcmp(dType,'curr') || strcmp(dType,'char') %attenuation correction for current/charge but not taus
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
        wt_Out(2:length(eachVel)+1,iCell) = num2cell(Iact);
    end
end
wt_Name = [{'wt', 'wt_Dist'}; thisName];
wt_Out = wt_Out(:,~cellfun(@isempty, wt_Out(1,:))); % clear out empty waves (where dist didn't match)
wt_Out = wt_Out(~cellfun(@isempty, wt_Out(:,1)),:);
wt_Name = wt_Name(~cellfun(@isempty, wt_Name(:,1)),:);

wt_Out = [[{'rampVel_wt'};num2cell(eachVel)] wt_Out]; % append stepSize wave


thisName = cell(0);
whichMRCs = fatMRCs;
fat_Out = cell(length(eachVel)+1,size(whichMRCs,1));

allVels = vertcat(fatMRCs{:,3});
allVels = allVels(:,1);
[~,~,eachVel] = sortRowsTol(allVels,tol,1);
if whichStim == 2
    eachVel = -eachVel;
end
if extFiltCorr
    eachVel(eachVel>19820)=19820;
    eachVel = unique(eachVel);
end

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist; %distance
        fat_Out{1,iCell} = [dType '_' thisName{iCell,1}];
        if extFiltCorr
            whichCorr=thisCell(:,1)>20000;
            thisCell(whichCorr,1) = 19820;
        end

        Iact = [];
        
        for iSize = 1:length(eachVel)
            stepSize = eachVel(iSize);
            whichStep = abs(roundVel(thisCell(:,1))-stepSize)<tol;
            if any(whichStep)
                Iact(iSize,1) = thisCell(whichStep,peakCol);
            else
                Iact(iSize,1) = nan;
            end
        end
        fat_Out(2:length(eachVel)+1,iCell) = num2cell(Iact);
    end
end

fat_Name = [{'fat', 'fat_Dist'}; thisName];
fat_Out = fat_Out(:,~cellfun(@isempty, fat_Out(1,:)));
fat_Out = fat_Out(~cellfun(@isempty, fat_Out(:,1)),:);
fat_Name = fat_Name(~cellfun(@isempty, fat_Name(:,1)),:);

fat_Out = [[{'rampVel_fat'};num2cell(eachVel)] fat_Out]; % append stepSize wave


xlswrite(fname,wt_Out,['wt_' dType]);
xlswrite(fname,fat_Out,['fat_' dType]);
xlswrite(fname,wt_Name,'wtStats');
xlswrite(fname,fat_Name,'fatStats');


clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName


% Normalize based on current/charge at given step size and write to xls
if normFlag && strcmp(dType,'curr') || strcmp(dType,'char')
    %WT
    which_Out = wt_Out;
    norm_Out = cell(size(which_Out));
    
    currMat = cell2mat(which_Out(2:end,2:end));
    whichRow = find(cellfun(@(x) abs(x-normVel)<tol, which_Out(2:end,1)));
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
    whichRow = find(cellfun(@(x) abs(x-normVel)<tol, which_Out(2:end,1)));
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

%% Calculate ratios
% re-save ant_Out as antOff/antOnCurr for the corresponding cases
on = antOnCurr;
off = antOffCurr;
out = on(1,:);

onVel = [on{2:end,1}];
offVel = -[off{2:end,1}];


for i = 1:length(onVel)
    out{i+1,1} = onVel(i);
    out(i+1,2:end) = cellfun(@(x,y) x./y, off(i+1,2:end), on(i+1,2:end),'un',0);
end

out = out(~cellfun(@isempty,out(:,1)),:);
out(1,:) = cellfun(@(x) regexprep(x,'stim1','ratio'),out(1,:),'un',0);

xlswrite(fname,out,[dType '_Ratio']);


%% Renormalize based on fit

% The point of the remaining code is to grab the data from Igor or from an Excel
% sheet after fitting individual I-d curves to get max/delta/xhalf from the
% Boltzmann fit, then normalize each I-d curve to its predicted max, and
% export back into an Igor-readable format for plotting.

% From AttCorrectedID_Curve.m

% From Igor, export the entire Igor file as a csv
%% Import fits from Igor table

fitStats = ImportMetaData(); %AttCorrectedI-D_fitstats.xlsx
fitHeaders = fitStats(1,:);
fitStats = fitStats(2:end,:);

%% Get names and match up stats
antStats = fitStats(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'Ant')), fitHeaders,'un',0)));
postStats = fitStats(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'Post')), fitHeaders,'un',0)));
antMaxIdx = find(~cellfun(@isempty,(regexp(fitHeaders(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'Ant')), fitHeaders,'un',0))),'[mM]ax'))));
postMaxIdx = find(~cellfun(@isempty,(regexp(fitHeaders(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'Post')), fitHeaders,'un',0))),'[mM]ax'))));
antStats = antStats(cellfun(@(x) ~isnumeric(x), antStats(:,1)),:);
postStats = postStats(cellfun(@(x) ~isnumeric(x), postStats(:,1)),:);

% this gives five columns: name, delta, distance, max, xhalf

% straight up taking the fits and recoloring them here won't work because 
% Igor fits are all 200 points long but the X wave
% is calculated, and is different length for each (i.e., some curves have
% X values from 0.5 to 11, some from 3 to 10, etc., and those ranges
% determine what the x values are for each fit separately).

% instead, I'm just going to do the normalization here and send the waves
% back to Igor, where I can average them like the other I-dCellFits

%%
eachVel = fitStats(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'stepSize')), fitHeaders,'un',0)));
eachVel = eachVel(cellfun(@(x) isnumeric(x) && ~isnan(x),eachVel));

whichSide = antStats;
thisCell = whichSide(:,1);
thisDataNorm = cell(0);

for i = 1:size(whichSide,1)
   thisMax = whichSide{i,antMaxIdx};
   whichData = cell2mat(cellfun(@(x) ~isempty(regexp(x,sprintf('^(?!fit).*%s',thisCell{i}))),fitHeaders,'un',0));
   thisData = fitStats(1:length(eachVel),whichData);
   thisData(cellfun(@isempty,thisData))={nan};
   thisDataNorm(:,i) = cellfun(@(x) x./thisMax, thisData,'un',0);
end

antNorm = thisDataNorm;
antHeaders = cellfun(@(x) sprintf('%s_Norm',x),thisCell,'un',0)';


% POSTERIOR
whichSide = postStats;
thisCell = whichSide(:,1);
thisDataNorm = cell(0);

for i = 1:size(whichSide,1)
   thisMax = whichSide{i,postMaxIdx};
   whichData = cell2mat(cellfun(@(x) ~isempty(regexp(x,sprintf('^(?!fit).*%s',thisCell{i}))),fitHeaders,'un',0));
   thisData = fitStats(1:length(eachVel),whichData);
   thisData(cellfun(@isempty,thisData))={nan};
   thisDataNorm(:,i) = cellfun(@(x) x./thisMax, thisData,'un',0);
end

postNorm = thisDataNorm;
postHeaders = cellfun(@(x) sprintf('%s_Norm',x),thisCell,'un',0)';

%%

xlswrite(fname,antHeaders,'antNorm');
xlswrite(fname,postHeaders,'postNorm');
xlswrite(fname,antNorm,'antNorm','A2');
xlswrite(fname,postNorm,'postNorm','A2');


