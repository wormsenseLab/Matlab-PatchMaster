% calcCalibs.m
%
% Allows you to select calibration step ranges and saves them to the given
% recording's calibration field in the ephysData struct.
%
%
%
%
%
% Created by Sammy Katta on 13 July 2016.

function ephysData = calcCalibs(ephysData)

% Import list of traces of ProbeCalib protocol (created by using
% ExcludeSweeps). e.g., CalibTraces104.xls
calibMetaData = ImportMetaData();
calibTracePicks = metaDataConvert(calibMetaData);
calibCells = calibMetaData(:,1);

% Read in photodiode traces (some early recordings also included current
% channel and thus had calibration signal in channel 3, not 2.
% If another channel is ever added, this will have to be redone.
for iCell = 1:length(calibCells)
    cellName = calibCells{iCell};
    
    % Load photodiode signal
    %TODO: Figure out if you actually need to save them all or if you can
    %just do this as plotData. (But you could output calibData/stepValues
    %if need be).
    if ephysData.(cellName).dataType{3} == 'V'
        calibData(:,iCell) = ephysData.(cellName).data{3,calibMetaData{iCell,2}};
    elseif ephysData.(cellName).dataType{2} == 'V'
        calibData(:,iCell) = ephysData.(cellName).data{2,calibMetaData{iCell,2}};
    else 
        disp('Photodiode calibration signal is neither channel 2 nor 3.');
    end

    
    clear handles;
    plotData = calibData(:,iCell);
    
    % Run GUI for selecting ranges (hit Begin Selection, then place data
    % cursors in pairs around your ranges. shift-Click for each new data
    % cursor, and click to move the most recently active cursor, or move
    % any cursor at a later time, but be sure to keep them paired with the
    % same start/end cursor).
    %TODO: Fix the datatip labels.
    handles = selectCalibSteps(plotData);
    
    % Pull out the mean value of the signal between the indices set by the
    % user's cursor placement. Store in stepValues with [recordings, steps]
    stepIdx = handles.cursorPoints(:,1);
    stepIdx = reshape(stepIdx,2,[])';
    for j = 1:size(stepIdx,1)
        stepValues(iCell,j) = mean(plotData(stepIdx(j,1):stepIdx(j,2)))';
    end
      
    close;

    ephysData.(cellName).calibration = [dVecAxial; stepValues(iCell,:)];


end



end