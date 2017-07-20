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

%Focus Adjust
calibDataPD{1} = ephysData.Calib04_syl_seal_focusadj.data{3,1};
calibDataPD{2} = ephysData.Calib04_syl_seal_focusadj.data{3,3};
calibDataPD{3} = ephysData.Calib04_syl_seal_focusadj.data{3,5};
calibDataPD{4} = ephysData.Calib04_syl_seal_focusadj.data{3,7};
calibDataPD{5} = ephysData.Calib04_syl_seal_focusadj.data{3,9};

calibDataProbe{1} = ephysData.Calib04_syl_seal_focusadj.data{2,2};
calibDataProbe{2} = ephysData.Calib04_syl_seal_focusadj.data{2,4};
calibDataProbe{3} = ephysData.Calib04_syl_seal_focusadj.data{2,6};
calibDataProbe{4} = ephysData.Calib04_syl_seal_focusadj.data{2,8};
calibDataProbe{5} = ephysData.Calib04_syl_seal_focusadj.data{2,10};


%% 
for iCell = 1:length(calibDataPD)
    clear handles;
    pdData = calibDataPD{iCell};
    probeData = calibDataProbe{iCell};
    
    % Run GUI for selecting ranges (hit Begin Selection, then place data
    % cursors in pairs around your ranges. shift-Click for each new data
    % cursor, and click to move the most recently active cursor, or move
    % any cursor at a later time, but be sure to keep them paired with the
    % same start/end cursor).
    %TODO: Fix the datatip labels.
    handlesPD = selectCalibSteps(pdData);
    handlesProbe = selectCalibSteps(probeData);
    
    % Pull out the mean value of the signal between the indices set by the
    % user's cursor placement. Store in stepValues with [recordings, steps]
    stepIdx = handlesPD.cursorPoints(:,1);
    stepIdx = reshape(stepIdx,2,[])';
    for j = 1:size(stepIdx,1)
        stepValuesPD(iCell,j) = mean(pdData(stepIdx(j,1):stepIdx(j,2)))';
    end
           
    close;
    
    stepIdx = handlesProbe.cursorPoints(:,1);
    stepIdx = reshape(stepIdx,2,[])';
    for j = 1:size(stepIdx,1)
        stepValuesProbe(iCell,j) = mean(probeData(stepIdx(j,1):stepIdx(j,2)))';
    end
    close;
    
end

% Ask where to save and write out the .xls file with calib mean values from
% PD signal channel
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save step values to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName,stepValuesPD,'PD');
xlswrite(fName,stepValuesProbe,'Probe');


% % Read in the xls file
% [filename, pathname] = uigetfile(...
%     {'*.xls;*.xlsx', 'Excel files';
%     '*.*', 'All files'}, ...
%     'Grab xls file with step values from calibrations');
% stepValuesPD = xlsread(fullfile(pathname,filename),'PD');
% stepValuesProbe = xlsread(fullfile(pathname,filename),'Probe');


probeSteps=stepValuesProbe;
pdSteps=stepValuesPD;
pdStepSizes = [0,10/3,20/3,10,40/3,50/3];
probeStepSizes = [0, 2.5, 5, 7.5, 10];
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


