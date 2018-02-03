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
whichChan = 2; %channel number for stimulus command
emptySeg = 4;
emptyFreq = 0;

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
        
        % Read from stimTree to find sine frequency and timepoints
        try thisTree = ephysData.(cellName).stimTree{thisSeries};
        catch      
            error('Please upload a stimTree for %s. (Run ImportPatchData with includeStim flag).\n',cellName);
        end
        
        chLoc = find(cellfun(@(x) ~isempty(x),thisTree(:,2)));
        
        try segLoc = chLoc(whichChan)+1:chLoc(whichChan+1)-1;
        catch
            segLoc = chLoc(whichChan)+1:size(thisTree,1);
        end
        
        thisFreq = 1/thisTree{chLoc(whichChan),2}.chSine_Cycle; % sine freq in Hz
        thisAmpl = thisTree{chLoc(whichChan),2}.chSine_Amplitude; % amplitude in commanded um
        
        for iSeg = 1:length(segLoc)
            segType (iSeg) = thisTree{segLoc(iSeg),3}.seClass; 
            segLength(iSeg) = thisTree{segLoc(iSeg),3}.seDuration * sf * 1e3; %length in samples
            
            if iSeg == 1
                segStart = [1];
            else
                segStart(iSeg) = segStart(iSeg-1)+segLength(iSeg-1);
            end
        end
        
        isSine = find(ismember(segType,3));
        if isempty(isSine)
            thisFreq = emptyFreq;
            thisAmpl = 0;
            isSine = emptySeg;
        end
        
        isSquare = find(ismember(segType,4));
        
        sineLoc = arrayfun(@(x) segStart(x):segStart(x)+segLength(x)-1, isSine, 'un', 0);
        %Check that this works as start/end points instead of all points
        squareLoc = arrayfun(@(x) [segStart(x), segStart(x)+segLength(x)-1], isSquare);

        for iSine = 1:length(sineLoc)
            sineTrace{iSine} = probeI(sineLoc{iSine});
        end
        
              
%NEXT: insert leak subtraction and pull sineTrace from leakSub instead
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


