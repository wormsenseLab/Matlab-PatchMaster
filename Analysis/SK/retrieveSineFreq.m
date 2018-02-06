% retrieveSineFreq.m
%
% This helper function retrieves the information about the frequency and
% amplitude of sine and square pulses for FrequencyAnalysis().
%
% sineParams: [sineStart sineEnd sineFreq sineAmpl]
% squareParams: [sqStart sqEnd sqAmpl nReps]

function [sineParams, squareParams] = retrieveSineFreq(ephysData, cellName, thisSeries, whichChan)

emptySeg = 4;
emptyFreq = 0;


try thisTree = ephysData.(cellName).stimTree{thisSeries};
catch
    return
end

chLoc = find(cellfun(@(x) ~isempty(x),thisTree(:,2)));

try segLoc = chLoc(whichChan)+1:chLoc(whichChan+1)-1;
catch
    segLoc = chLoc(whichChan)+1:size(thisTree,1);
end

sineFreq = 1/thisTree{chLoc(whichChan),2}.chSine_Cycle; % sine freq in Hz
sineAmpl = thisTree{chLoc(whichChan),2}.chSine_Amplitude; % amplitude in commanded um
sqAmpl = thisTree{chLoc(whichChan),2}.chSquare_PosAmpl;
sf = 1/thisTree{1,1}.stSampleInterval;

for iSeg = 1:length(segLoc)
    segType (iSeg) = thisTree{segLoc(iSeg),3}.seClass;
    segLength(iSeg) = thisTree{segLoc(iSeg),3}.seDuration * sf; %length in samples
    
    if iSeg == 1
        segStart = [1];
    else
        segStart(iSeg) = segStart(iSeg-1)+segLength(iSeg-1);
    end
end

isSine = find(ismember(segType,3));
if isempty(isSine)
    sineFreq = emptyFreq;
    sineAmpl = 0;
    isSine = emptySeg;
end

isSquare = find(ismember(segType,4));
sineParams = arrayfun(@(x) [segStart(x), segStart(x)+segLength(x)-1], isSine, 'un', 0);
sineParams = vertcat(sineParams{:});

sineParams(:,3) = repmat(sineFreq, size(sineParams,1),1);
sineParams(:,4) = repmat(sineAmpl, size(sineParams,1),1);

squareParams = arrayfun(@(x) [segStart(x), segStart(x)+segLength(x)-1], isSquare, 'un', 0);

squareParams = repmat(squareParams,2,1);  
squareParams = vertcat(squareParams{:});
squareParams(:,3) = repmat(sqAmpl, size(squareParams,1),1);
squareParams(:,4) = ones(size(squareParams,1),1);

% Separate on and off stimuli for use in findMRCs(). Assume because this is
% a square pulse that step command is instantaneous (i.e.,
% stepStart:stepEnd is one point).
for iSq = 1:2:size(squareParams,1)
    squareParams(iSq,2) = squareParams(iSq,1)+1;
end
for iSq = 2:2:size(squareParams,1)
    squareParams(iSq,1) = squareParams(iSq,2)+1;
    squareParams(iSq,2) = squareParams(iSq,1)+1;
    squareParams(iSq,3) = -squareParams(iSq,3);
end


end