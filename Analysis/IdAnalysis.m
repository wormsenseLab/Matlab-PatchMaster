% IdAnalysis.m
%
% Created by Sammy Katta on 20-May-2015.

function mechPeaks = IdAnalysis(ephysData, allCells)

% keyboard;

mechPeaks = cell(length(allCells),1);
sf = 5; %sampling frequency, in kHz
stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
baseLength = baseTime*sf;
smoothWindow = 5; % n timepoints for moving average window for findPeaks

% Load Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx;*.csv;*.txt', 'All spreadsheets';
    '*.xls;*.xlsx', 'Excel files';
    '*.csv', 'Comma-separated value files';
    '*.txt', 'Tab-delimited text files';
    '*.*', 'All files'}, ...
    'Pick spreadsheets with list of files/series to analyze',...
    'MultiSelect', 'on');
fName = fullfile(pathname, filename);
[~,~,mechTracePicks] = xlsread(fName);
for i = 1:length(mechTracePicks)
    mechTracePicks{i,4} = str2num([mechTracePicks{i,3}]);
end
mechTracePicks = mechTracePicks(:, [1 2 4]);


%     allSteps = [];
%     allOns = [];
%     allOffs = [];
%     allOnTaus = [];
%     allOffTaus = [];




% Find applicable series and check against list of included series/traces
% (this allows a cross-check on the protocol name) before analyzing
% Values for traces not on the list will be stored as NaN.
for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    protName = 'WC_Probe';
    allSeries = find(strcmp(protName,ephysData.(cellName).protocols));
    protName = 'WC_ProbeLarge';
    allSeries = [allSeries find(strcmp(protName,ephysData.(cellName).protocols))];
    protName = 'WC_ProbeSmall';
    allSeries = [allSeries find(strcmp(protName,ephysData.(cellName).protocols))];
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
    allSizes = [];
    allLeakSub = [];
    allStarts = [];
    allEnds = [];
    
    nSeries = length(allSeries);
    
    for iSeries = 1:nSeries
        
        % Carry out analysis if this series is on the list
        try pickedTraces = pickedSeries{[pickedSeries{:,1}]==allSeries(iSeries),2};
        catch
            continue
        end
        
        probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
        % convert command V to um, at 0.408 V/um
        stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
        nSteps = size(stimComI,2);
        
        % Initialize arrays as NaNs (can be placed in array). You can
        % later use nanmean to take the mean while ignoring NaN values.
        stepSize = NaN(nSteps,1);
        stepStarts = NaN(nSteps,1);
        stepEnds = NaN(nSteps,1);
        leakSubtract = NaN(length(stimComI),nSteps);
        %             onPeaks = NaN(nSteps,1);
        %             offPeaks = NaN(nSteps,1);
        %             onTau = NaN(nSteps,1);
        %             offTau = NaN(nSteps,1);
        %
        % Find timepoint in sweep at which probe starts pushing (on) and when it
        % goes back to neutral position (off). Use that timepoint to find nearby
        % peak current
        
        for iStep = 1:nSteps
            
            % Analyze only desired traces within this series
            if any(pickedTraces == iStep)
                
                % Figure out timepoints when stimulus command for step 
                % starts and ends
                stepOn = stimComI(:,iStep) - mean(stimComI(1:10*sf,iStep));
                stepStart = find(stepOn > stepThresh);
                stepLength = length(stepStart);
                if stepLength == 0
                    continue
                end
                stepStarts(iStep) = stepStart(1);
                stepEnds(iStep) = stepStart(1) + stepLength;
                stepSize(iStep) = mean(stepOn(stepStarts(iStep)+1:stepEnds(iStep)-1));
                              
                % Subtract leak/baseline
                leak = mean(probeI(1:100,iStep));
                leakSubtract(:,iStep) = probeI(:,iStep) - leak;
                
            end
        end
        
        % Concatenate to the complete list of step sizes and
        % leak-subtracted traces across series for this recording
        allSizes = [allSizes; stepSize];
        allStarts = [allStarts; stepStarts];
        allEnds = [allEnds; stepEnds];
        allLeakSub = [allLeakSub; leakSubtract'];
        clear leakSubtract;
              
    end
    
    % Round step size to nearest 0.1 to allow grouping data by step size
    allSizes = round(allSizes*10)/10;
    
    % Sort by size and take start/end indices of the data for each size
    [sortedSizes, sortIdx] = sort(allSizes);
    [eachSize,sizeStartIdx,~] = unique(sortedSizes,'first');
    [~,sizeEndIdx,~] = unique(sortedSizes,'last');
    
    sortedStarts = allStarts(sortIdx);
    sortedEnds = allEnds(sortIdx);
    sortedLeakSub = allLeakSub(sortIdx,:);
    
    % TODO: Store nReps as endIdx-StartIdx for each step size
    
    % Use start index for the start and end times, assuming they don't
    % change within a given step size (or whatever grouping you are using;
    % should work with different step lengths/intervals as well).
    startsBySize = sortedStarts(sizeStartIdx);
    endsBySize = sortedEnds(sizeStartIdx);
    
    % Use start and end indices for each step size to take the mean of the
    % leak-subtracted trace corresponding to that step size. Then smooth
    % and find peaks near the step times.
    for iSize = 1:sum(~isnan(eachSize))
        meansBySize(iSize,:) = mean(sortedLeakSub(sizeStartIdx(iSize):sizeEndIdx(iSize),:));
        
        % Smooth data with a moving average for peak finding
        smooMean = -smooth(meansBySize(iSize,:)',smoothWindow,'moving');
        % Set threshold based on noise of the first 100ms of the trace
        % (i.e., size of signal needed to be seen above that noise)
        pkThresh(iSize) = 1.5*thselect(smooMean(1:100*sf),'rigrsure');
        
        % Find MRC peaks if they exist at the onset of the step, otherwise
        % set peak amplitude as NaN
        % TODO: Consider whether you want to just take pkLoc and find the
        % amplitude from the original trace, with a baseline subtracted
        % from the immediate pre-stimulus time period (original plan).
        [pk, pkLoc] = findpeaks(smooMean(startsBySize(iSize)-100:startsBySize(iSize)+300),...
            'minpeakheight',pkThresh(iSize));
        if ~isempty(pk)
            pkOn(iSize) = max(pk);
            pkLoc = pkLoc(pk==max(pk));
            pkOnLoc(iSize) = pkLoc(1) + startsBySize(iSize)-100; %account for start position
        else
            pkOn(iSize) = 0;
            pkOnLoc(iSize) = nan;
        end
        
        % Find MRC peaks at the offset of the step
        [pk, pkLoc] = findpeaks(smooMean(endsBySize(iSize)-100:endsBySize(iSize)+300),...
            'minpeakheight',pkThresh(iSize));
        if ~isempty(pk)
            pkOff(iSize) = max(pk);
            pkLoc = pkLoc(pk==max(pk));
            pkOffLoc(iSize) = pkLoc(1) + endsBySize(iSize)-100;
        else
            pkOff(iSize) = 0;
            pkOffLoc(iSize) = nan;
        end
        
        
    end
    
 
    % Get baseline for each step by grabbing the mean of the 150 points before
    % the probe displacement.
    
    baseProbeI_on = mean(sortedLeakSub(:,startsBySize-baseLength:startsBySize),2);
    baseProbeI_off = mean(sortedLeakSub(:,endsBySize-baseLength:endsBySize),2);
    
    
    % Find the peak current for the on step and the off step for this sweep
    
% onSubtract(stepStart+1:stepStart+500) == max(onSubtract(stepStart+1:stepStart+500)));
    peakOnLoc = peakOnLoc(1) + stepStart;
    
    peakOffLoc = find(offSubtract(stepEnd+1:stepEnd+500) == max(offSubtract(stepEnd+2:stepEnd+500)));
    peakOffLoc = peakOffLoc(1)+stepEnd;
    
    onPeaks(iStep) = -onSubtract(peakOnLoc)*1E12;
    offPeaks(iStep) = -offSubtract(peakOffLoc)*1E12;
    
    [~,onFitInd] = min(abs(onSubtract(peakOnLoc:75*sf+peakOnLoc)-(onSubtract(peakOnLoc)/(2*exp(1)))));
    [~,offFitInd] = min(abs(offSubtract(peakOffLoc:75*sf+peakOffLoc)-(offSubtract(peakOffLoc)/(2*exp(1)))));
    
    onFitTime = onFitInd/sf; % seconds
    onT = 0:1/sf:onFitTime;
    offFitTime = offFitInd/sf; % seconds
    offT = 0:1/sf:offFitTime;
    
    
    onFit = fit(onT',onSubtract(peakOnLoc:peakOnLoc+onFitInd),'exp1');
    offFit = fit(offT',onSubtract(peakOffLoc:peakOffLoc+offFitInd),'exp1');
    %     plot(capFit,t,ICt(intStart:intStart+minInd));
    
    % Calculate time constant in seconds for calculation
    onTau(iStep) = -1/onFit.b;
    offTau(iStep) = -1/offFit.b;
    % TODO: Output taus
    % TODO: Fit with an alpha function instead of exp1 to get tau1 and tau2
    
    
    
    
    
    
    
    mechPeaks{iCell} = [allSteps allOns allOffs allOnTaus allOffTaus];
    
    % TODO: Figure out how to fit this to the four-parameter sigmoidal
    % function used in O'Hagan: @(X,a,b,c,d) ((a-d)/(1+((X/c)^b)))+d
    % Using optimtool? fmincon? nlinfit if you add the statistics toolbox.
    
    
    
    
    
end

end