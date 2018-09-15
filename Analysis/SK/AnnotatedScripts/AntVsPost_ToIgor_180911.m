% Tau vs distance curves for anterior vs. posterior

% Running all the way through allows user to pull out tau data from steps, 
% export to Igor, fit a Boltzmann there, re-import and normalize and export
% back to Igor for final plotting and fitting. 

% Unlike AntVsPost_VoltageAttenuation, this does not correct the current or
% charge for voltage errors based on distance, and was used for non-current
% data.

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
stimPosition = {'posterior'};

wormPrep = {'dissected'};
cellDist = [40 90]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = 2.5; % frequency of low-pass filter for stimulus command
includeFlag = 1; 

posteriorDistCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', includeFlag);

stimPosition = {'anterior'};

anteriorDistCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included',includeFlag);


clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep excludeCells;

% This results in a file with three columns: cell ID, series, and sweeps 
% that have been approved for analysis use.

%% Visual/manual exclusion of bad sweeps 
% This is mainly based on excluding sweeps with leak > 10pA, but also
% sweeps where the recording was lost partway through or some unexpected
% source of noise was clearly at play.

protList ={'WC_Probe';'WC_ProbeSmall';'WC_ProbeLarge';'NoPrePulse'};
matchType = 'full';

ExcludeSweeps(ephysData, protList, anteriorDistCells, 'matchType', matchType);
ExcludeSweeps(ephysData, protList, posteriorDistCells, 'matchType', matchType);

%% Find MRCs 
% antAllSteps.xlsx and postAllSteps.xlsx from 180810
protList ={'WC_Probe';'NoPre'};
matchType = 'first';

sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

anteriorMRCs = IdAnalysis(ephysData,protList,anteriorDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1);

posteriorMRCs = IdAnalysis(ephysData,protList,posteriorDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1);

clear protList sortSweeps matchType

%% Pull out and combine relevant data for plot

% Voltage attenuation data comes from the length constant fitting done in
% Igor, which gives a voltage attenuation factor at the location of the
% stimulus site for each recording, based on the calculated length
% constant.

% Grab the current at the selected step size (e.g., 10um) for each
% recording, grab the cell-stimulus distance for that recording, and match
% the recording name against the voltage attenuation table to find the
% attenuation factor for that recording.
whichMRCs = anteriorMRCs;
thisAtt = attenuationData(:,[2 8 10]);
distCol = 12;
peakCol = 8;     % 8 for tau act, 9 for tau decay
distVPeak = [];
stepSize = 10;

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    whichStep = round(thisCell(:,1)) == stepSize;
    if any(whichStep)
        thisName = whichMRCs{iCell,1};
        hasAtt = strcmp(thisName,thisAtt(:,1));
        
        if any(hasAtt) && thisAtt{hasAtt,2}
            distVPeak(iCell,:) = [thisCell(whichStep,[distCol peakCol]) thisAtt{hasAtt,3}];
        else
            distVPeak(iCell,:) = [thisCell(whichStep,[distCol peakCol]) nan];
        end
        
    end
    
end

distVPeak_Ant = distVPeak;

whichMRCs = posteriorMRCs;
distVPeak = [];


for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    whichStep = round(thisCell(:,1)) == stepSize;
    if any(whichStep)
        thisName = whichMRCs{iCell,1};
        hasAtt = strcmp(thisName,thisAtt(:,1));
        
        if any(hasAtt)
            distVPeak(iCell,:) = [thisCell(whichStep,[distCol peakCol]) thisAtt{hasAtt,2}];
        else
            distVPeak(iCell,:) = [thisCell(whichStep,[distCol peakCol]) nan];
        end
        
    end
    
end

distVPeak_Post = distVPeak;

clear a iCell thisCell whichMRCs whichStep thisName hasAtt distVPeak

%% Export for Igor fitting of Boltzmann to each recording

% This is the same as in AntVPost_VoltageAttenuation_180803 except it
% doesn't apply any attenuation correction. Useful for pulling out the taus

eachSize = [0.5 1 1.5 3 4 5 6 7 8 9 10 11 12]';
distLimits = [40 90];
% limit to same average distance for anterior and posterior

whichMRCs = anteriorMRCs;
thisAtt = attenuationData(:,[2 8 10]);
thisName = cell(0);
ant_Out = cell(length(eachSize)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        ant_Out{1,iCell} = thisName{iCell,1};
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
    thisCell = whichMRCs{iCell,3};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist; %distance
        post_Out{1,iCell} = thisName{iCell,1};
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

% Set the filename
fname = 'PatchData/taus_antVpost_distLim(180911).xls';

xlswrite(fname,ant_Out,'antAct');
xlswrite(fname,post_Out,'postAct');
xlswrite(fname,ant_Name,'antActStats');
xlswrite(fname,post_Name,'postActStats');


clear iCell thisCell Iact whichMRCs whichStep hasAtt thisAtt Vc Ena thisName


%% Renormalize

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