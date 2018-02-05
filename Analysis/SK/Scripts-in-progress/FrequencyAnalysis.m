% FrequencyAnalysis.m

function output = FrequencyAnalysis(ephysData, ephysMetaData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('ephysMetaData', @(x) iscell(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));

p.parse(ephysData, ephysMetaData, protList, varargin{:});

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

sinePeaks = cell(length(allCells),1);

baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline for leak subtraction
steadyTime = 100; % length of time (ms) to use for calculating steady state mean at end of sine
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
whichChan = 2; %channel number for stimulus command

for iCell = 1:length(allCells)
    
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
        
        
        % Subtract the baseline leak current for each series
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtractCell = num2cell(leakSubtract',2);
                

        clear thisFreq;
        % Find location, frequency, and amplitude of sine, and locations of
        % flanking square pulses
        [thisFreq, thisAmpl, sineLoc, squareLoc] = ...
            retrieveSineFreq(ephysData, cellName, thisSeries, whichChan);
        
        if ~exist(thisFreq,'var')
            fprintf('Please upload a stimTree for %s. (Run ImportPatchData with includeStim flag).\n',cellName);
        end

        % For later: when using multiple sines in a given sweep        
        % for iSine = 1:length(sineLoc)
        %     sineTrace{iSine} = leakSubtract(sineLoc{iSine});
        % end
        
        sineTrace = leakSubtract(sineLoc{1});
        
        steadyStateI = mean(sineTrace(end-steadyTime*sf:end));
        
        
%NEXT: take average of beginning and end of sine section (x ms? x% of
%trace? x cycles?) and save to output cell array
%output: sine freq, amplitude, duration, peak at start, beginning average, end
%average, tau of fit of smoothed trace
%NEXT: findpeaks on square pulses (on and off based on squareLoc) and
%output: sine freq, amplitude, on1, on2, off1, off2, interval between
%squares, interval between sweeps

%THEN: sort stim profiles based on freq/ampl/dur, and take averages per
%frequency.

%THEN: output freq dependence and adaptation like mechPeaks, with cell
%name, trace, cell array containing each recording's average output
        

    end
end

end


