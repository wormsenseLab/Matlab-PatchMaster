%Define data to work with

% calibDataPD{1} = ephysData.E1D.data{2,5};
% calibDataPD{2} = ephysData.E1F.data{2,3};
% calibDataPD{3} = ephysData.E1G.data{2,5};

% calibDataProbe{1} = ephysData.E1D.data{3,6};
% calibDataProbe{2} = ephysData.E1F.data{3,2};
% calibDataProbe{3} = ephysData.E1G.data{3,4};

% %Model Cell
% calibDataPD{1} = ephysData.Calib01_mod_WC.data{3,4};
% calibDataPD{2} = ephysData.Calib02_mod_seal.data{3,2};
% calibDataPD{3} = ephysData.Calib03_syl_seal.data{3,7};
% calibDataPD{4} = ephysData.Calib03_syl_seal.data{3,10};
% calibDataPD{5} = ephysData.Calib05_syl_seal_pip.data{3,1};
% 
% calibDataProbe{1} = ephysData.Calib01_mod_WC.data{2,5};
% calibDataProbe{2} = ephysData.Calib02_mod_seal.data{2,3};
% calibDataProbe{3} = ephysData.Calib03_syl_seal.data{2,8};
% calibDataProbe{4} = ephysData.Calib03_syl_seal.data{2,11};
% calibDataProbe{5} = ephysData.Calib05_syl_seal_pip.data{3,3};

% %Focus Adjust
% calibDataPD{1} = ephysData.Calib04_syl_seal_focusadj.data{3,1};
% calibDataPD{2} = ephysData.Calib04_syl_seal_focusadj.data{3,3};
% calibDataPD{3} = ephysData.Calib04_syl_seal_focusadj.data{3,5};
% calibDataPD{4} = ephysData.Calib04_syl_seal_focusadj.data{3,7};
% calibDataPD{5} = ephysData.Calib04_syl_seal_focusadj.data{3,9};
% 
% calibDataProbe{1} = ephysData.Calib04_syl_seal_focusadj.data{2,2};
% calibDataProbe{2} = ephysData.Calib04_syl_seal_focusadj.data{2,4};
% calibDataProbe{3} = ephysData.Calib04_syl_seal_focusadj.data{2,6};
% calibDataProbe{4} = ephysData.Calib04_syl_seal_focusadj.data{2,8};
% calibDataProbe{5} = ephysData.Calib04_syl_seal_focusadj.data{2,10};

% allCells = fieldnames(ephysData);
% allCells = allCells(~cellfun(@isempty,allCells));


calibPDTracePicks = ImportMetaData(); 
calibPDTracePicks = metaDataConvert(calibPDTracePicks);
calibPD2TracePicks = ImportMetaData(); % AllWCStepsTo104TracePicks.xls
calibPD2TracePicks = metaDataConvert(calibPD2TracePicks);
calibProbeTracePicks = ImportMetaData(); % AllWCStepsTo104TracePicks.xls
calibProbeTracePicks = metaDataConvert(calibProbeTracePicks);
calibPDTracePicks= [calibPDTracePicks;calibPD2TracePicks];
clear calibPD2TracePicks

calibPDCells = calibPDTracePicks(:,1);
calibProbeCells = calibProbeTracePicks(:,1);

for iCell = 1:13
    calibDataPD{iCell} = ephysData.(calibPDCells{iCell}).data{2,calibPDTracePicks{iCell,2}};
    calibDataProbe{iCell} = ephysData.(calibProbeCells{iCell}).data{3,calibProbeTracePicks{iCell,2}};
end

for iCell = 14:length(calibPDCells)
    calibDataPD{iCell} = ephysData.(calibPDCells{iCell}).data{3,calibPDTracePicks{iCell,2}};
    calibDataProbe{iCell} = ephysData.(calibProbeCells{iCell}).data{3,calibProbeTracePicks{iCell,2}};
end

%% Set timepoints for averages

%NEXT:Run this and select regions, also on wormsubtract

for iCell = 38:length(calibDataPD)
    clear handles;
    pdData3 = calibDataPD{iCell};
    probeData3 = calibDataProbe{iCell};
    
    % Run GUI for selecting ranges (hit Begin Selection, then place data
    % cursors in pairs around your ranges. shift-Click for each new data
    % cursor, and click to move the most recently active cursor, or move
    % any cursor at a later time, but be sure to keep them paired with the
    % same start/end cursor).
    %TODO: Fix the datatip labels.
    handlesPD = selectCalibSteps(pdData3,calibPDCells{iCell});
    handlesProbe = selectCalibSteps(probeData3,calibProbeCells{iCell});
    
    % Pull out the mean value of the signal between the indices set by the
    % user's cursor placement. Store in stepValues with [recordings, steps]
    stepIdx = handlesPD.cursorPoints(:,1);
    stepIdx = reshape(stepIdx,2,[])';
    for j = 1:size(stepIdx,1)
        stepValuesPD(iCell,j) = mean(pdData3(stepIdx(j,1):stepIdx(j,2)))';
    end
           
    close;
    
    stepIdx = handlesProbe.cursorPoints(:,1);
    stepIdx = reshape(stepIdx,2,[])';
    for j = 1:size(stepIdx,1)
        stepValuesProbe(iCell,j) = mean(probeData3(stepIdx(j,1):stepIdx(j,2)))';
    end
    close;

%     seriesStimuli = newStepFind(1, probeData3, 10, 10);
    
end

clear pdData3 handlesPD handlesProbe stepIdx iCell j

% Ask where to save and write out the .xls file with calib mean values from
% PD signal channel
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save step values to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName,stepValuesPD,'PD');
xlswrite(fName,stepValuesProbe,'Probe');

%% Find commanded z for probe steps
clear cmdStepValues stepValuesCmd

for iCell = 1:length(calibProbeTracePicks)
    calibDataProbeStep = ephysData.(calibProbeCells{iCell}).data{2,calibProbeTracePicks{iCell,2}};   
    sf = ephysData.(calibProbeCells{iCell}).samplingFreq{calibProbeTracePicks{iCell,2}}/1000;
    cmdStepValues{iCell,1} = newStepFind(1,calibDataProbeStep,sf);
    for j = 1:size (cmdStepValues{iCell},1)
        stepValuesCmd(iCell,j) = cmdStepValues{iCell}(j,4);
    end
end

%Pad with zero at beginning bc stepfinding doesn't include baseline atm
stepValuesCmd = [zeros(size(cmdStepValues,1),1) stepValuesCmd];

%% Load in commanded z for PD steps
%Create excel sheet with probe step protocols since this data is not saved
%in the dat file, for each series being analyzed (like tracePicks but w
%step sizes in the last column), and read it in/vectorize it.

calibPDSteps = ImportMetaData(); % CalibTracesPD_Steps.xls
calibPDSteps = metaDataConvert(calibPDSteps);
calibPD2Steps = ImportMetaData(); % CalibTracesPD2_Steps.xls
calibPD2Steps = metaDataConvert(calibPD2Steps);
calibPDSteps= [calibPDSteps;calibPD2Steps];
clear calibPD2Steps

%% Load in image tracking data from csv and extract for plotting

calibBase = ImportMetaData(); % Poker\PDCalibrationDatabase.xlsx

ephysData = ImportCalibImgTracks(ephysData, calibBase);


for iCell = 1:length(calibProbeTracePicks)
    try calibImgTrace = ephysData.(calibProbeCells{iCell}).imgTrack(:,1);   
    catch
        disp(sprintf('%s was not used.', calibProbeCells{iCell}));
        continue
    end
    [~,imgStepLocs]=findpeaks(abs(diff(calibImgTrace)),'MinPeakProm',0.05);
    imgStepMeans = arrayfun(@(x) mean(calibImgTrace(x-10:x)),imgStepLocs);
    imgStepMeans(end+1) = mean(calibImgTrace(end-10:end));
    imgStepValues{iCell,1} = imgStepMeans';
    stepValuesImg(iCell,:) = imgStepMeans';
end

%NEXT: also output imgStepValues as double (put in try catch to avoid
%future issues) and create plots vs. pz-z-cmd and pz-PD-signal

clear imgStepMeans iCell calibImgTrace
%% Plots 

%Subtract offset
stepValuesProbeZ = stepValuesProbe - repmat(stepValuesProbe(:,1),1,11);
stepValuesPDZ = stepValuesPD - repmat(stepValuesPD(:,1),1,11);

%Plot stepValuesCmd (Vcmd) vs. imgStepValues (to show if linear)
%Add second x-axis to this plot that scales Vcmd --> z-cmd-pz
%TODO: Add mean trace
% stepValuesImg = ImportMetaData();
% stepValuesImg = reshape([stepValuesImg{:}],4,11);

a = stepValuesCmd(~stepValuesImg(:,1)==0,:);
b = stepValuesImg(~stepValuesImg(:,1)==0,:);

plot(stepValuesCmd(~stepValuesImg(:,1)==0,:)'/0.408,...
    stepValuesImg(~stepValuesImg(:,1)==0,:)')
hold on
xlim([0 15])
fig.ax(1) = gca; % current axes
ax1_pos = fig.ax(1).Position; % position of first axes
fig.ax(2) = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XLim', [0 15*0.408]);
fig.ax(2).XColor = 'r';
set(fig.ax(2),'YTick',[])
set(fig.ax(2),'YColor','k')
xlabel('Probe command (V)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(1));
hold on
plot([0,15],[0,15],'k--')
ylabel('Imaged horiz movement (um)'); %based on PD signal
xlabel('Probe command (um)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(2));


%Plot this new z-cmd-PD vector vs. stepValuesPD (V-PD)
%DECIDE just the last four or all or some subset? Include load?
fig.f1 = figure(); 
% fig.f2 = figure();
for i=15:length(calibPDSteps)
    if ~isempty(strfind(calibPDSteps{i,1},'unload')) || ~isempty(strfind(calibPDSteps{i,1},'same'))
        figure(fig.f1);
        plot(calibPDSteps{i,3},stepValuesPDZ(i,1:length(calibPDSteps{i,3})),'r');
    else
        figure(fig.f1);
        plot(calibPDSteps{i,3},stepValuesPDZ(i,1:length(calibPDSteps{i,3})),'b');
    end
    hold on
end
xlabel('Commanded PD movement (mm)')
ylabel('PD signal (V)')

%Plot stepValuesCmd(scaled to z-cmd-pz) vs. stepValuesProbe (V-pz)

fig.f3 = figure(); 
% fig.f4 = figure();
for i=15:length(calibPDSteps)
    if ~isempty(strfind(calibPDSteps{i,1},'unload')) || ~isempty(strfind(calibPDSteps{i,1},'same'))
        figure(fig.f3);
   plot(stepValuesCmd(i,:)/0.408,stepValuesProbeZ(i,:),'r')
    else
        figure(fig.f3);
   plot(stepValuesCmd(i,:)/0.408,stepValuesProbeZ(i,:),'b')
    end
    hold on
end
hold on
xlim([0 15])
fig.ax(1) = gca; % current axes
ax1_pos = fig.ax(1).Position; % position of first axes
fig.ax(2) = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XLim', [0 15*0.408]);
fig.ax(2).XColor = 'r';
set(fig.ax(2),'YTick',[])
set(fig.ax(2),'YColor','k')
xlabel('Probe command (V)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(1));
ylabel('PD signal (V)'); %based on PD signal
xlabel('Probe command (um)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(2));

%Calculate/interpolate z-act-pz from z-cmd-PD vs. V-PD given V-pz
%Plot stepValuesCmd(scaled to z-cmd-pz) vs. z-act-pz

for i = 15:length(calibPDSteps)
stepInterpProbe(i,:) = interp1(stepValuesPDZ(i,:),...
    calibPDSteps{i,3}*1000/60,...
    stepValuesProbeZ(i,:),'linear','extrap');
end

fig.f5=figure();
for i=15:length(calibPDSteps)
    if ~isempty(strfind(calibPDSteps{i,1},'unload')) || ~isempty(strfind(calibPDSteps{i,1},'same'))
        plot(stepValuesCmd(i,:)'/0.408,stepInterpProbe(i,:)','r')
    else
        plot(stepValuesCmd(i,:)'/0.408,stepInterpProbe(i,:)','b')
    end
    hold on
end
xlim([0 15])
fig.ax(1) = gca; % current axes
ax1_pos = fig.ax(1).Position; % position of first axes
fig.ax(2) = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XLim', [0 15*0.408]);
fig.ax(2).XColor = 'r';
set(fig.ax(2),'YTick',[])
set(fig.ax(2),'YColor','k')
xlabel('Probe command (V)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(1));
hold on
plot([0,15],[0,15],'k--')
ylabel('Calculated horiz movement (um)'); %based on PD signal
xlabel('Probe command (um)'); %based on geometry (assuming 17deg angle)
axes(fig.ax(2));


%Plot imgStepValues vs. z-act-pz


% Bring in new unloaded vs. loaded data

%% Look at centered vs not centered

centerCells = {'pzvspd_load_A_center', 'pzvspd_loadB_center', 'pzvspd_loadC_center', 'pzvspd_loadD_center'}';
edgeCells = {'pzvspd_load_A3', 'pzvspd_load_B', 'pzvspd_loadC', 'pzvspd_loadD'}';

fig.f6 = figure();
for i = 1:length(centerCells)
    plot();
    
end


%%

% % Read in the xls file
% [filename, pathname] = uigetfile(...
%     {'*.xls;*.xlsx', 'Excel files';
%     '*.*', 'All files'}, ...
%     'Grab xls file with step values from calibrations');
% stepValuesPD = xlsread(fullfile(pathname,filename),'PD');
% stepValuesProbe = xlsread(fullfile(pathname,filename),'Probe');


probeSteps=stepValuesProbe;
pdSteps=stepValuesPD;
% pdStepSizes = [0,10/3,20/3,10,40/3,50/3];
% probeStepSizes = [0, 2.5, 5, 7.5, 10];
% probeStepSizes = 0.143:0.143:1;

for i = 1:5
pdInterp(i,:) = interp1(pdStepSizes,pdSteps(i,:),probeStepSizes,'linear');
end

probeStepsZ = probeSteps - repmat(probeSteps(:,1),1,5);
pdStepsZ = pdSteps - repmat(pdSteps(:,1),1,6);

for i = 1:5
pdInterpZ(i,:) = interp1(pdStepSizes,pdStepsZ(i,:),probeStepSizes,'linear');
end

for i = 1:5
pdInterpZProbe(i,:) = interp1(pdStepsZ(i,:),pdStepSizes,probeStepsZ(i,:),'linear');
end

figure()
plot(pdInterpZ',probeStepsZ','b')
plotfixer
hold on
plot([-4,2],[-4, 2],'k--')
xlabel('Interpolated PD movement signal (V)');ylabel('Probe movement signal (V)');

figure()
plot(probeStepSizes,pdInterpZProbe') 
plotfixer
hold on;
plot([0,15],[0,15],'k--')
ylabel('Calculated horiz movement (um)'); %based on PD signal
xlabel('Commanded horiz probe movement (um)'); %based on geometry (assuming 17deg angle)


