% I-d and Q-d curves for velocity steps (corrected for space clamp)
% Current or charge vs. displacement, and vs. distance 

% This version uses TrapRate protocols, which have both the on and the off
% at the same velocity (except for one particular velocity where
% Patchmaster refused to do the off at the same velocity), in a trapezoid
% shape.

% At the moment, this is limited to anterior recordings, filtered at both
% 2.5 or 5 kHz (because we can see the actual velocity).

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

wormPrep = {'dissected'};
cellDist = [40 200]; % stimulus/cell distance in um
resistCutoff = '<250'; % Rs < 250 MOhm
extFilterFreq = [2.5 5]; % frequency of low-pass filter for stimulus command
includeFlag = 1; 

stimPosition = {'anterior'};

velocityCells = FilterRecordings(ephysData, ephysMetaData,...
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

protList ={'TrapRate'};
matchType = 'first';

ExcludeSweeps(ephysData, protList, velocityCells, 'matchType', matchType);

%% Find MRCs 
% antAllSteps.xlsx and postAllSteps.xlsx from 180810
protList ={'TrapRate'};
matchType = 'first';

sortSweeps = {'velocity','velocity','magnitude','magnitude'};

velocityMRCs = IdAnalysis(ephysData,protList,velocityCells,'num','matchType',matchType, ...
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
whichMRCs = velocityMRCs;
thisAtt = attenuationData(:,[2 8 10]);
distCol = 12;
peakCol = 6; % 6 for peak current, 11 for integrated current/charge
distVPeak = [];
stepVel = 11360;

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    whichVel = round(thisCell(:,1)) == stepVel;
    if any(whichVel)
        thisName = whichMRCs{iCell,1};
        hasAtt = strcmp(thisName,thisAtt(:,1));
        
        if any(hasAtt) && thisAtt{hasAtt,2}
            distVPeak(iCell,:) = [thisCell(whichVel,[distCol peakCol]) thisAtt{hasAtt,3}];
        else
            distVPeak(iCell,:) = [thisCell(whichVel,[distCol peakCol]) nan];
        end
        
    end
    
end

distVPeak_Ant = distVPeak;

clear a iCell thisCell whichMRCs whichStep thisName hasAtt distVPeak

%% Make correction to I for space clamp error
% (Note: this is an approximation, assuming that the majority of the
% current is happening at the stimulus site, which we know is not entirely
% true, and channels at different points along the neurite will experience
% different voltages.

% Now calculate what the actual voltage was at the stimulus site based on
% the voltage attenuation factor. In this case, command voltage was -60mV
% for all steps. 

% attenuation factor = Vm'/Vc so Vm' = Vc * attenuation factor
 
% Then calculate what current would've been based on sodium reversal
% potential. 
% For IC2 solution used here, Ena = +94mV.
 
% Im' = Im * ((Vc-Ena)/(Vm'-Ena))
 
Vc = -0.06; %in V
Ena = 0.094; % in V
Im = distVPeak_Ant(:,2);
Vatt = distVPeak_Ant(:,3);

Vm = Vc * Vatt;
distVPeak_Ant(:,4) = (Im * (Vc-Ena)) ./ (Vm-Ena);

%% Plot anterior and posterior on separate plots

figure(); axh(1) = subplot(2,1,1);
scatter(distVPeak_Ant(:,1),distVPeak_Ant(:,4));
xlabel('Distance from cell body (um)');
ylabel(sprintf('Current @%dum (pA)',stepVel))
text(100,200,'Anterior');

for i = 1:length(axh)
    xlim(axh(i),[0 200]);
    ylim(axh(i),[0 200]);
end

%% Correct all sizes and export for Igor fitting of Boltzmann to each recording

% Set the filename
fname = 'PatchData/attCorrectedVelCurves(180913).xls';

allVel = vertcat(velocityMRCs{:,3});
eachVel = uniquetol(allVel(:,1),12,'DataScale',1);
% the value of 12 is empirically chosen as being larger than the largest
% difference between velocities that could be binned (usually 10 after roundVel) and smaller than the
% difference between velocities that shouldn't (i.e., 64 vs. 80).
distLimits = [40 200];
% limit to same average distance for anterior and posterior


Vc = -0.06; %in V
Ena = 0.094; % in V

whichMRCs = velocityMRCs;
thisAtt = attenuationData(:,[2 8 10]);
thisName = cell(0);
vel_Out = cell(length(eachVel)+1,size(whichMRCs,1));

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
    if thisDist <= distLimits(2) && thisDist >= distLimits(1)
        thisName{iCell,1} = whichMRCs{iCell,1}; %name
        thisName{iCell,2} = thisDist;
        vel_Out{1,iCell} = thisName{iCell,1};
        hasAtt = strcmp(thisName{iCell,1},thisAtt(:,1));
        Iact = [];
        
        for iVel = 1:length(eachVel)
            stepVel = eachVel(iVel);
            whichVel = abs(thisCell(:,1) - stepVel) < 12; %find closest velocity
            if any(whichVel)
                
                if any(hasAtt) && thisAtt{hasAtt,2} %if attenuation calc exists and not omitCell
                    Vm = Vc * thisAtt{hasAtt,3};
                    Im = thisCell(whichVel,peakCol);
                    Iact(iVel,1) = (Im * (Vc-Ena)) ./ (Vm-Ena);
                else
                    Iact(iVel,1) = nan;
                end
                
            else
                Iact(iVel,1) = nan;
            end
        end
        vel_Out(2:length(eachVel)+1,iCell) = num2cell(Iact);
    end
end
vel_Name = [{'Anterior', 'Ant_Dist'}; thisName];
vel_Out = vel_Out(:,~cellfun(@isempty, vel_Out(1,:))); % clear out empty waves (where dist didn't match)
vel_Name = vel_Name(~cellfun(@isempty, vel_Name(:,1)),:);

vel_Out = [[{'velocity'};num2cell(eachVel)] vel_Out]; % append stepSize wave

% thisName = cell(0);
% whichMRCs = posteriorMRCs;
% post_Out = cell(length(eachSize)+1,size(whichMRCs,1));
% 
% for iCell = 1:size(whichMRCs,1)
%     thisCell = whichMRCs{iCell,3};
%     thisDist = mean(thisCell(:,distCol)); % check if cell distance is in range
%     if thisDist <= distLimits(2) && thisDist >= distLimits(1)
%         thisName{iCell,1} = whichMRCs{iCell,1}; %name
%         thisName{iCell,2} = thisDist; %distance
%         post_Out{1,iCell} = thisName{iCell,1};
%         Iact = [];
%         
%         for iSize = 1:length(eachSize)
%             stepVel = eachSize(iSize);
%             whichStep = round(thisCell(:,1)*2)/2 == stepVel; %round to nearest 0.5
%             if any(whichStep)
%                 Iact(iSize,1) = thisCell(whichStep,peakCol);
%             else
%                 Iact(iSize,1) = nan;
%             end
%         end
%         post_Out(2:length(eachSize)+1,iCell) = num2cell(Iact);
%     end
% end
% 
% post_Name = [{'Posterior', 'Post_Dist'}; thisName];
% post_Out = post_Out(:,~cellfun(@isempty, post_Out(1,:)));
% post_Name = post_Name(~cellfun(@isempty, post_Name(:,1)),:);


xlswrite(fname,vel_Out,'vel');
xlswrite(fname,vel_Name,'velStats');

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