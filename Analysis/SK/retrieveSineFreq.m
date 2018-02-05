% retrieveSineFreq.m
%
% This helper function retrieves the information about the frequency and
% amplitude of sine and square pulses for FrequencyAnalysis().
%

function [thisFreq, thisAmpl, sineLoc, squareLoc] = retrieveSineFreq(ephysData, cellName, thisSeries, whichChan)

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

sineLoc = arrayfun(@(x) [segStart(x), segStart(x)+segLength(x)-1], isSine, 'un', 0);
%Check that this works as start/end points instead of all points
squareLoc = arrayfun(@(x) [segStart(x), segStart(x)+segLength(x)-1], isSquare, 'un', 0);

end