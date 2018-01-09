% NonStatNoiseAnalysis.m
% 
% 
% 
% Created by Sammy Katta on 5th Jan 2018.

function [noiseTraces] = NonStatNoiseAnalysis(ephysData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));


p.parse(ephysData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

% Allow for pre-filtered subsets of "allCells" (from FilterRecordings) to 
% be used - otherwise take names from imported metadata sheet.
if isempty(allCells)
    allCells = unique(mechTracePicks(:,1));
end

averagingWindow = 8; %sliding window size (n sweeps) for noise analysis

stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um



for iCell = 1:length(allCells)
    cellName = allCells(iCell)
   
    
    allStim = cell(0);
    allLeakSub = cell(0);
    
    % Double check that the list of series given matches the indices of the
    % protocol names specified in the input.
    
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,protList,'MatchType',matchType);
    
    nSeries = length(allSeries);
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
    if nSeries == 0 || isempty(pickedSeries)
        continue
    end
    
    for iSeries = 1:nSeries
        thisSeries = allSeries(iSeries);
        
        % Carry out analysis if this series is on the list
        try pickedTraces = pickedSeries{[pickedSeries{:,1}]==thisSeries,2};
        catch
            continue % if it's not on the list, go on to next series in for loop
        end
        
        probeI = ephysData.(cellName).data{1,thisSeries}(:,pickedTraces);
        stimComI = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces); %in V, not um
        % sampling frequency in kHz
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        dataType = ephysData.(cellName).dataunit{1,thisSeries};
        nSweeps = size(stimComI,2);

                
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtractCell = num2cell(leakSubtract',2);
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);

        %NEXT: for iSweep=1:averagingWindow/2:nSweeps-1, take average of x
        %leak-subtracted sweeps within the window, then subtract average
        %from each individual sweep.
        %NEXT: use seriesStimuli location/sweep number to set time
        %boundaries for where to look at the response
        %NEXT: save variance (subtracted trace) and mean current for each
        %window, plot vs. each other for all sweeps
        
        windowMeans = cell(0);
        meanSubtract = cell(0);
        
        for iSweep = 1:averagingWindow/2:nSweeps-averagingWindow+1
            try theseSweeps = probeI(:,iSweep:iSweep+averagingWindow-1);
                windowMeans{iSweep} = mean(theseSweeps,2);
            catch
                continue
            end
            
            meanSubtract{iSweep} = theseSweeps-repmat(windowMeans{iSweep},1,averagingWindow);
            
        end
        
        
    end
    
end

end