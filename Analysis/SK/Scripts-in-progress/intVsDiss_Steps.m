% IntactAnalysis.m 
% I-d and Q-d curves for anterior dissected vs. anterior intact
% preparations

%%Import data and metadata
ephysData = ImportPatchData('incl',1);
projects = {'FAT';'SYM'};
ephysData = FilterProjectData(ephysData, projects);
clear projects;

%fix a couple messed up names
ephysData.FAT029 = ephysData.FAT029s;
ephysData.FAT017 = ephysData.FAT017e001;
ephysData.FAT164 = ephysData.FAT164001;
ephysData = rmfield(ephysData,{'FAT017e001';'FAT029s';'FAT164001'});

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
cellDist = [40 90]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = [2.5 5]; % frequency of low-pass filter for stimulus command
includeFlag = 1; 

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
ExcludeSweeps(ephysData, protList, intactDistCells, 'matchType', matchType);

%% Find MRCs 
% antAllSteps.xlsx and postAllSteps.xlsx from 180810
protList ={'WC_Probe';'NoPre'};
matchType = 'first';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

[anteriorMRCs, antStim] = IdAnalysis(ephysData,protList,dissectedDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

posteriorMRCs = IdAnalysis(ephysData,protList,intactDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1,'subZeroCharge',1);

clear protList sortSweeps matchType



%% Correct all sizes and export for Igor fitting of Boltzmann to each recording

% Set the filename
fname = 'PatchData/antVPost_off_allDist_subQ(190215).xls';


eachSize = [0.5 1 1.5 3 4 5 6 7 8 9 10 11 12]';
distLimits = [40 200]; % limit to same average distance for anterior and posterior
distCol = 12;
whichStim = 2; %on or off
noCorr = 0; % 1 to skip space clamp voltage attenuation correction
Vc = -0.06; %in V
Ena = 0.094; % in V


whichMRCs = anteriorMRCs;
dType = 'char'; 

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

thisName = cell(0);
whichMRCs = posteriorMRCs;
post_Out = cell(length(eachSize)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,whichStim+2};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist; %distance
        post_Out{1,iCell} = [dType '_' thisName{iCell,1}];
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
        post_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
    end
end

post_Name = [{'Posterior', 'Post_Dist'}; thisName];
post_Out = post_Out(:,~cellfun(@isempty, post_Out(1,:)));
post_Name = post_Name(~cellfun(@isempty, post_Name(:,1)),:);

post_Out = [[{'stepSize'};num2cell(eachSize)] post_Out]; % append stepSize wave



xlswrite(fname,ant_Out,['ant_' dType]);
xlswrite(fname,post_Out,['post_' dType]);
xlswrite(fname,ant_Name,'antStats');
xlswrite(fname,post_Name,'postStats');


clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName


%% Normalize based on current/charge at given step size and write to xls
% Normalizes within each condition

stepSize = -10; % must be negative for off

%Anterior
which_Out = ant_Out;
norm_Out = cell(size(which_Out));

currMat = cell2mat(which_Out(2:end,2:end));
whichRow = find(cellfun(@(x) x==stepSize, which_Out(2:end,1)));
normRow = currMat(whichRow,:);
normMat = bsxfun(@rdivide,currMat,normRow);
norm_Out(2:end,2:end) = num2cell(normMat);
norm_Out(:,1) = which_Out(:,1);
norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);


ant_Norm = norm_Out;


% Posterior
which_Out = post_Out;
norm_Out = cell(size(which_Out));

currMat = cell2mat(which_Out(2:end,2:end));
whichRow = find(cellfun(@(x) x==stepSize, which_Out(2:end,1)));
normRow = currMat(whichRow,:);
normMat = bsxfun(@rdivide,currMat,normRow);
norm_Out(2:end,2:end) = num2cell(normMat);
norm_Out(:,1) = which_Out(:,1);
norm_Out(1,2:end) = cellfun(@(x) ['norm_' x],which_Out(1,2:end),'un',0);

post_Norm = norm_Out;

xlswrite(fname,ant_Norm,['ant_' dType '_Norm']);
xlswrite(fname,post_Norm,['post_' dType '_Norm']);

clear currMat which Row normRow normMat norm_Out

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

xlswrite(fname,out,['ant_' dType '_Ratio']);


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
eachSize = fitStats(:,cell2mat(cellfun(@(x) ~isempty(regexp(x,'stepSize')), fitHeaders,'un',0)));
eachSize = eachSize(cellfun(@(x) isnumeric(x) && ~isnan(x),eachSize));

whichSide = antStats;
thisCell = whichSide(:,1);
thisDataNorm = cell(0);

for i = 1:size(whichSide,1)
   thisMax = whichSide{i,antMaxIdx};
   whichData = cell2mat(cellfun(@(x) ~isempty(regexp(x,sprintf('^(?!fit).*%s',thisCell{i}))),fitHeaders,'un',0));
   thisData = fitStats(1:length(eachSize),whichData);
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
   thisData = fitStats(1:length(eachSize),whichData);
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


%% Plot representative traces for spatial dynamics figure

% check distance and nReps for each recording
antDists = cell2mat(cellfun(@(x) mean(x(:,12:13)),anteriorMRCs(:,3),'un',0));
postDists = cell2mat(cellfun(@(x) mean(x(:,12:13)),posteriorMRCs(:,3),'un',0));

figure();
plot(posteriorMRCs{10,2}');
title('FAT224 (posterior)');
ylim([-8e-11 2e-11]);
cmapline('ax',gca,'colormap','copper');

figure();
plot(anteriorMRCs{3,2}');
title('FAT143 (anterior)');
ylim([-8e-11 2e-11]);
cmapline('ax',gca,'colormap','copper');

plotfixer

%% Separate and write representative traces for channel-sim-distance figure

% anterior trace FAT 105 at 1, 5, 10um mean traces separated by on/off
% stimulus, from [stimStart-250, stimEnd+250]ms.
onStim = 751;
offStim = 2251;
sf = 5; %kHz
boxTime = 150*sf; %ms
tVec = (-150*sf:150*sf)'/sf/1000; %s

reps(:,1:3) = anteriorMRCs{4,2}([2 5 9],onStim-boxTime:onStim+boxTime)';
reps(:,4:6) = anteriorMRCs{4,2}([2 5 9],offStim-boxTime:offStim+boxTime)';


