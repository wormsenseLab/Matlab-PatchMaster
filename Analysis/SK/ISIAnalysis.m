% ISIAnalysis.m
% 
% 

function ISIPeaks = ISIAnalysis(ephysData, allCells, protName)

% keyboard;
sf = 5; %sampling frequency, in kHz
stepThresh = 0.05; % step detection threshold in um, could be smaller
ISIPeaks = cell(length(allCells),1);
smoothWindow = 5;

for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    allSeries = find(strcmp(protName,ephysData.(cellName).protocols));
    
    allOns = [];
    allOffs = [];
    
    nSeries = length(allSeries);
    tStart = zeros(1,nSeries);
    
    for iSeries = 1:nSeries
        probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
        % convert command V to um, at 0.408 V/um for current setup(5-15-15)
        stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
        tStart(1,iSeries) = ephysData.(cellName).startTimes{allSeries(iSeries)};
        
        nSteps = size(stimComI,2);
        
        stepSize = NaN(nSteps,1);
        leakSubtract = NaN(nSteps,length(stimComI));

        pkThresh = zeros(nSteps,1);
        pkOn = zeros(nSteps,1);
        pkOff = zeros(nSteps,1);
        pkOnLoc = NaN(nSteps,1);
        pkOffLoc = NaN(nSteps,1);
        onsetTau = NaN(nSteps,1);
        offsetTau = NaN(nSteps,1);

        for iStep = 1:nSteps
            % Figure out points when stimulus command for step starts and ends
            stepOn = stimComI(:,iStep) - mean(stimComI(1:10*sf,iStep));
            stepStart = find(stepOn > stepThresh); 
            stepLength = length(stepStart);
            if stepLength == 0
                continue
            end
            stepStart = stepStart(1);
            stepEnd = stepStart + stepLength;
            stepSize(iStep) = mean(stepOn(stepStart+1:stepEnd-1));
            
            % Get baseline for each step by grabbing the mean of the first
            % 100 points, then subtract that leak from the trace.
            leak = mean(probeI(1:100,iStep));
            leakSubtract(iStep,:) = probeI(:,iStep) - leak;
            
%             % Find the peak current for the on step and the off step for this sweep
%             % Two options for ignoring the stimulus artifact (which shows up
%             % even with unpatched pipette):
%             % Start looking at stepStart +2 (+1 too early, +3 hits actual peak)
%             % Swap subtract order and don't use abs(onSubtract) for finding max
% 
%             
%             onSubtract = baseProbeI_on - probeI(:,iStep);
%             
%             
%             peakOnLoc = find(onSubtract(stepStart+1:stepStart+500) == max(onSubtract(stepStart+1:stepStart+500)));
%             peakOnLoc = peakOnLoc(1) + stepStart;
%             
%             offSubtract = baseProbeI_off - probeI(:,iStep);
%             peakOffLoc = find(offSubtract(stepEnd+1:stepEnd+500) == max(offSubtract(stepEnd+2:stepEnd+500)));
%             peakOffLoc = peakOffLoc(1)+stepEnd;
%             
%             onPeaks(iStep) = -onSubtract(peakOnLoc);
%             offPeaks(iStep) = -offSubtract(peakOffLoc);
            
            
            
                %%%%%%%%
                
        %TODO: Put the peak finding into a separate function and pass in
        %the info for stepStart vs stepEnd.
        
            % Smooth data with a moving average for peak finding
        smooMean = -smooth(leakSubtract(iStep,:)',smoothWindow,'moving');
        % Set threshold based on noise of the first 100ms of the trace
        % (i.e., size of signal needed to be seen above that noise)
        pkThresh(iStep) = 1.5*thselect(smooMean(1:100*sf),'rigrsure');
        
        % Find MRC peaks if they exist at the onset of the step, otherwise
        % set peak amplitude as NaN. Calculate decay constant tau based on
        % single exponent fit for onset and offset currents.
        % TODO: Consider whether you want to just take pkLoc and find the
        % amplitude from the original trace, with a baseline subtracted
        % from the immediate pre-stimulus time period (original plan).
        [pk, pkLoc] = findpeaks(smooMean(stepStart-100:stepStart+100),...
            'minpeakheight',pkThresh(iStep));
        if ~isempty(pk)
            
            pkOn(iStep) = max(pk)*1E12;
            pkLoc = pkLoc(pk==max(pk));
            pkOnLoc(iStep) = pkLoc(1) + stepStart-100; %account for start position
            
            % Find time for current to decay to 2/e of the peak or 75ms
            % after the peak, whichever comes first. Use that for fitting
            % the single exponential. Fit the unsmoothed mean trace.
            [~,onFitInd] = min(abs(leakSubtract(iStep,pkOnLoc(iStep):75*sf+pkOnLoc(iStep))...
                - (leakSubtract(pkOnLoc(iStep))/(2*exp(1)))));
            
            onFitTime = onFitInd/sf; % seconds
            onT = 0:1/sf:onFitTime;
            
            onFit = fit(onT',leakSubtract(iStep,pkOnLoc(iStep):pkOnLoc(iStep)+onFitInd)','exp1');
            onsetTau(iStep) = -1/onFit.b;
            
        end
        
        % Find MRC peaks at the offset of the step
        [pk, pkLoc] = findpeaks(smooMean(stepEnd-100:stepEnd+100),...
            'minpeakheight',pkThresh(iStep));
        if ~isempty(pk)
            pkOff(iStep) = max(pk)*1E12;
            pkLoc = pkLoc(pk==max(pk));
            pkOffLoc(iStep) = pkLoc(1) + stepEnd-100;
            
%             [~,offFitInd] = min(abs(leakSubtract(iStep,pkOffLoc(iStep):75*sf+pkOffLoc(iStep))...
%                 - (leakSubtract(pkOffLoc(iStep))/(2*exp(1)))));
%             
%             offFitTime = offFitInd/sf; % seconds
%             offT = 0:1/sf:offFitTime;
%             
%             offFit = fit(offT',leakSubtract(iStep,pkOffLoc(iStep):pkOffLoc(iStep)+offFitInd)','exp1');
%             offsetTau(iStep) = -1/offFit.b;
% 
        end
        
    
    
    %%%%%%%%%%%%%%%%
            
        end
        
        allOns = [allOns pkOn];
        allOffs = [allOffs pkOff];

        
        
        
    end
    

    
    ISIPeaks{1,iCell} = allOns;
    ISIPeaks{2,iCell} = allOffs;
    ISIPeaks{3,iCell} = tStart;
   
end

% Adjust start times so beginning of first ISI series is time 0.
ISIPeaks(3,:) = cellfun(@(x) x-x(1),ISIPeaks(3,:),'UniformOutput',0);

% On and off adjusted to zero
ISIPeaks(4,:) = cellfun(@(x) x-x(1),ISIPeaks(1,:),'UniformOutput',0);
ISIPeaks(5,:) = cellfun(@(x) x-x(1),ISIPeaks(2,:),'UniformOutput',0);

end