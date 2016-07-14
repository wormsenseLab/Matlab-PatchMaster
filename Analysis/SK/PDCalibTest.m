% PDCalibTest.m
%
%

%TODO: Read in list of calibrations? Or autosearch for ProbeCalib protocol
%and perform this on all instances of that?
% nCalibs = 3;

% Temporarily defining calibration protocol locations, later, auto-detect
% and associate each with the given ephysData recording rather than just
% the stepValues variable.
% calibData(:,1) = ephysData.TestCalibrationWorm1.data{3,4};
% calibData(:,2) = ephysData.Worm_B.data{3,11};
% calibData(:,3) = ephysData.Worm_C.data{3,9};
% calibData(:,1) = ephysData.PutUnloadCalib.data{3,44};
% % calibData(:,2) = ephysData.PutUnloadCalib.data{3,46};
% calibData(:,2) = ephysData.Calibrated.data{2,11};
% calibData(:,3) = ephysData.Calibrated.data{2,12};

calibMetaData = ImportMetaData(); %CalibTraces104.xls
calibTracePicks = metaDataConvert(calibMetaData);
calibCells = calibMetaData(:,1);

for iCell = 1:4
    calibData(:,iCell) = ephysData.(calibCells{iCell}).data{3,calibMetaData{iCell,2}};
end

for iCell = 5:length(calibCells)
    calibData(:,iCell) = ephysData.(calibCells{iCell}).data{2,calibMetaData{iCell,2}};
    
end

%% Select ranges from calibration protocol
for iCell = 1:length(calibCells)
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
    
end

%% Split up wt and mutants
genotype = cell(length(calibCells),2);
for iCell=1:length(calibCells)
    stepValues(iCell,:) = stepValues(iCell,:)-stepValues(iCell,1);
    
    genotype(iCell,1) = calibCells(iCell);
    genotype(iCell,2) = ephysRecordingBase(strcmp(ephysRecordingBase(:,1),calibCells(iCell)),2);    
end
wtCalibCells = calibCells(strcmp(genotype(:,2),'TU2769'));
fatCalibCells = calibCells(strcmp(genotype(:,2),'GN381'));

wtStepValues = stepValues(ismember(calibCells,wtCalibCells),:);
fatStepValues = stepValues(ismember(calibCells,fatCalibCells),:);


%% Calculate calibration curve and convert photodiode readings into actual step sizes

%TODO: Have dVec included as input in calibTracePicks (metaDataConvert has
%been modified to allow it as column 4).
dVec = [0 2.5 5 7.5 10];
% dcmap = [0.4 0.4 0.4; 0.7 0.7 0.7; 1 0 0];
% set(gca, 'ColorOrder', dcmap)
% hold on
dVecAxial = dVec / cos(deg2rad(17));
plot(dVecAxial,wtStepValues(3:end,:),'k')
hold on
plot(dVecAxial,fatStepValues(1:end,:),'r')
plotfixer;

baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline

for iCell = 1:length(calibCells)
   allPDSteady = []; 
   allTraceIDs = []; 
   allSeriesIDs = [];
   
   cellName = calibCells{iCell};
   calibSteps(:,:,iCell) = ephysData.(cellName).data{1,calibMetaData{iCell,2}}; 
   
   allSeries = matchProts(ephysData,cellName,...
       {'WC_Probe','Probe_CC'},'MatchType','full');
   nSeries = length(allSeries);
   pickedSeries = calibTracePicks(find(strcmp(cellName,calibTracePicks(:,1))),[2,3]);
   
   for iSeries = 1:nSeries
       thisSeries = allSeries(iSeries);
       
       % Carry out analysis if this series is on the list
       try pickedTraces = pickedSeries{[pickedSeries{:,1}]==thisSeries,2};
       catch
           continue % if it's not on the list, go on to next series in for loop
       end
       
%        probeI = ephysData.(cellName).data{1,thisSeries};
       % convert command V to um, at 0.408 V/um
%TODO: include stimComI and call findSteps and use allSizes to sort rather
%than traceID, because that won't work if you include Small/Large protocols
       photodiodeV = ephysData.(cellName).data{3,thisSeries};
       % sampling frequency in kHz
       sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
%        dataType = ephysData.(cellName).dataunit{1,thisSeries};
       nSteps = size(photodiodeV,2);
       
       pdChange = SubtractLeak(photodiodeV, sf, 'BaseLength', baseTime);
       allPDSteady = [allPDSteady pdChange(1750:2249,:)];
       allTraceIDs = [allTraceIDs pickedTraces];
       allSeriesIDs = [allSeriesIDs repmat(thisSeries, length(pickedTraces),1)];
   end
   
   [sortedTraceID, sortIdx] = sort(allTraceIDs); % sort by trace number, assuming size is same
   [eachSize,sizeStartIdx,~] = unique(sortedTraceID,'first');
   [~,sizeEndIdx,~] = unique(sortedTraceID,'last');
   nSizes = sum(~isnan(eachSize));
   sortedPDSteady = allPDSteady(:,sortIdx);
   sortedPDAverage = mean(sortedPDSteady,1);
   sortedSeriesID = allSeriesIDs(sortIdx); % combine with sortedTraceID to ID traces used for mean
   
   meansBySize = NaN(nSizes,length(sortedPDSteady));

   for iSize = 1:nSizes
       if sizeEndIdx(iSize)-sizeStartIdx(iSize)>0
           meansBySize(iSize,:) = mean(sortedPDSteady(:,sizeStartIdx(iSize):sizeEndIdx(iSize)),2);
       else
           meansBySize(iSize,:) = sortedPDSteady(:,sizeStartIdx(iSize):sizeEndIdx(iSize));
       end 
   end
   
   if ~nSeries ==0
       pdAverage(:,iCell) = mean(meansBySize,2);
   else
       pdAverage(:,iCell) = nan;
   end
       


end

%% Interpolate 

wtPDAverages = pdAverage(:,ismember(calibCells,wtCalibCells));
fatPDAverages = pdAverage(:,ismember(calibCells,fatCalibCells));

for iCell = 3:length(wtPDAverages)
try wtStepCurves(iCell,:) = interp1(-wtStepValues(iCell,:), dVecAxial, -wtPDAverages(:,iCell)', 'linear','extrap');
catch
    wtStepCurves(iCell,:) = 0;
end

end

for iCell = 1:length(fatPDAverages)
try fatStepCurves(iCell,:) = interp1(-fatStepValues(iCell,:), dVecAxial, -fatPDAverages(:,iCell)', 'linear','extrap');
catch
    fatStepCurves(iCell,:) = 0;
end

end

%% Plot actual step curves
plot([1 3 5 7 9 11], wtStepCurves,'b');
hold on
plot([1 3 5 7 9 11], fatStepCurves,'g');
plot([1 3 5 7 9 11], [1 3 5 7 9 11],':k')
plotfixer;
% Pull out wt and fat PD step data from WC_Probe or Probe_CC steps (just
% stick with the 1,3,5,7,9,11 for now)
% wt: 82,83,84,86,96,97
% fat: 89,90,92,94,102,104
% StepCalibs79to104.xls
% 1750:2249 for 5kHz sample
% For a given recording, put all WC/CC steps in array [length 6 serieses]
% Take mean across serieses, then mean across 1750:2249 of length
% Should give 6 values per recordingn
% (Also take mean across length first and plot across series as check)
% Use those 6 values to interp actual displacement, plot vs commanded

