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

baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um

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
        

%NEXT: read in and match prots from ExcludeSweeps output as usual
%For each group, search stimTree for sine frequency and determine which
%segment has the sine/the timepoints of the sine.
%Q: how to include the 0_Sine control? Use similar stim profile sorting and
%keep times constant?
%Calculate the average of the last n ms of the sine window

%Optional analysis: look for flanking square pulses and calculate on/on and
%off/off ratio.


    end
end

end


