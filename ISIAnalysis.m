% ISIAnalysis.m
% 
% 

function ISIPeaks = ISIAnalysis(ephysData, allCells)

% keyboard;
sf = 5; %sampling frequency, in kHz
stepThresh = 0.05; % step detection threshold in um, could be smaller
ISIPeaks = cell(length(allCells),1);

for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    protName = 'WC_Probe5';
    allSeries = find(strcmp(protName,ephysData.(cellName).protocols));
    
    allOns = [];
    allOffs = [];
    
    nSeries = length(allSeries);
    
    for iSeries = 1:nSeries
        probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
        % convert command V to um, at 0.408 V/um for current setup(5-15-15)
        stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
        
        nSteps = size(stimComI,2);
        
        stepSize = zeros(nSteps,1);
        onPeaks = zeros(nSteps,1);
        offPeaks = zeros(nSteps,1);

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
            
            % Get baseline for each step by grabbing the mean of the 250 points before
            % the probe displacement.
            baseProbeI_on = mean(probeI(stepStart-150+1:stepStart-1,iStep));
            baseProbeI_off = mean(probeI(stepEnd-150+1:stepEnd-1,iStep));

            
            % Find the peak current for the on step and the off step for this sweep
            % Two options for ignoring the stimulus artifact (which shows up
            % even with unpatched pipette):
            % Start looking at stepStart +2 (+1 too early, +3 hits actual peak)
            % Swap subtract order and don't use abs(onSubtract) for finding max
            onSubtract = baseProbeI_on - probeI(:,iStep);
            peakOnLoc = find(onSubtract(stepStart+1:stepStart+500) == max(onSubtract(stepStart+1:stepStart+500)));
            peakOnLoc = peakOnLoc(1) + stepStart;
            
            offSubtract = baseProbeI_off - probeI(:,iStep);
            peakOffLoc = find(offSubtract(stepEnd+1:stepEnd+500) == max(offSubtract(stepEnd+2:stepEnd+500)));
            peakOffLoc = peakOffLoc(1)+stepEnd;
            
            onPeaks(iStep) = -onSubtract(peakOnLoc);
            offPeaks(iStep) = -offSubtract(peakOffLoc);
        end
        
        allOns = [allOns; onPeaks];
        allOffs = [allOffs; offPeaks];

    end
    
    ISIPeaks{iCell} = [allOns allOffs];
   
end