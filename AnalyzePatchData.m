% AnalyzePatchData.m
% 
% TODO: Figure out import for RecordingDatabase.xlsx, add columns for OC,WC
% series numbers of interest and pull those.
% TODO: Pass in series numbers from spreadsheet for analysis, rather than
% analyzing all series with that protocol name (but check protocol name of
% those series and pass error if not).

%% Import Data

% Don't forget to run sigTOOL first!
ephysData = ImportPatchData();

% Keep only data with given project prefixes/names.
projects = {'FAT'};

ephysData = FilterProjectData(ephysData, projects);

clear projects
%% Analyze capacity transient for C, Rs, and tau

ephysData = CtAnalysis(ephysData);

%% Temporary: Assign numbers of IVq series to look at

% TODO: Read these in from a separate file.

allCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT027';'FAT028';'FAT029';'FAT030';'FAT031'; 'FAT032'};


% List which sets of ct_ivqs to use for on cell (row 1)/whole cell (row 2)
% calculations for the above selection of cells. (i.e., if the second  
% ct_ivq protocol run for that recording was the one you want to use for  
% whole-cell, put "2" in row 1). Make sure it has three sequential ivq
% pgfs.
protStart = [1 1 1 4 1 1 1 1 1 1 1; ...
             4 4 4 7 4 7 6 4 4 4 4];
         
% Which Rs value to use for IV Rs correction
protRs = ceil(protStart(2,:)./3);
protRs(7) = 3;
         
for i = 1:length(allCells)
    ephysData.(allCells{i}).protOC = protStart(1,i);
    ephysData.(allCells{i}).protWC = protStart(2,i);    
    ephysData.(allCells{i}).protRs = protRs(i);
end
clear i protStart allCells protRs
%% Process voltage steps

% allCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT027';'FAT028';'FAT029';'FAT030';'FAT031';'FAT032'};
% 
% allIVs = IVAnalysis(ephysData,allCells);

wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031'};
fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032'};

testingSplit = IVAnalysis(ephysData,wtCells,fatCells);
wtIVs = testingSplit{1};
fatIVs = testingSplit{2};

clear testingSplit wtCells fatCells

%% Correct voltage steps based on series resistance and plot as scatter

wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031'};
fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032'};

wtIVs = IVRsCorrection(ephysData,wtIVs,wtCells);
fatIVs = IVRsCorrection(ephysData,fatIVs,fatCells);

wtI = [];
wtV = [];
fatI = [];
fatV = [];

for i = 1:length(wtCells)
    cellName = wtCells{i};
    wtV = [wtV;wtIVs.(cellName).actualV'];
    wtI = [wtI;wtIVs.(cellName).meanI'];
end

for i = 1:length(wtCells)
    cellName = fatCells{i};
    fatV = [fatV;fatIVs.(cellName).actualV'];
    fatI = [fatI;fatIVs.(cellName).meanI'];
end

figure()
hold on;
scatter(wtV*1E3,wtI*1E12);
scatter(fatV*1E3,fatI*1E12,'d');
plotfixer;
% FAT027 and FAT030 are the weird ones

%% Plot mechanically evoked currents in response to single steps
% IdAnalysis

wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031'};
fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032'};

mechPeaksWT = IdAnalysis(ephysData,wtCells);
mechPeaksFat = IdAnalysis(ephysData,fatCells);

%% Plot single MRC sets
% Draw stim protocol for MRCs
dt = 0.2; % ms, based on sampling frequency (5kHz in current WC_Probe)
tVec = 0:dt:750-dt;
dVec = zeros(length(tVec),12);

% for i = 1:6
%     dVec(50/dt:250/dt,i)=1+2*(i-1);
% end
% 
% figure()
% plot(tVec,dVec,'r');
% plotfixer;

% Get data
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,1);
% plot(tVec,toPlot/1E-12);
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,2);
% plot(tVec,toPlot/1E-12);
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,3);
% plot(tVec,toPlot/1E-12);
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,4);
% plot(tVec,toPlot/1E-12);
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,5);
% plot(tVec,toPlot/1E-12);
% figure()
% toPlot = ephysData.FAT030.data{1,59}(:,6);
% plot(tVec,toPlot/1E-12);

figure()
toPlot = ephysData.FAT032.data{1,18};
plot(tVec,toPlot/1E-12,'b');
plotfixer;

%% ISI validation

% Validate that interstimulus interval is long enough for full recovery
% from stimulus by running a step with the ISI of choice 16-32x. Check to
% see whether on/off currents amplitude is steady over time.
allCells = {'FAT032'};

ISIPeaksTest = ISIAnalysis(ephysData,allCells);

figure();
plot(ISIPeaksTest{1}*1E12);
plotfixer;

%% Dt 
allCells = {'FAT032'};

OnDtTest = OnDtAnalysis(ephysData,allCells);

figure()
scatter(OnDtTest{1}(:,1),OnDtTest{1}(:,3))
hold on;
scatter(OnDtTest{1}(:,1),OnDtTest{1}(:,2),'d')
% scatter(OnDtTest{1}(:,1),OnDtTest{1}(:,4),'s')


%% Look at current clamp
allCells = fieldnames(ephysData);

for iCell = 1:length(allCells)
    cellName = allCells{iCell}; %split into project name and cell numbers when feeding input
    
    % UPDATE: after June/July 2014 (FAT025), Patchmaster has separate pgfs for
    % 'OC_ct_neg' and 'WC_ct_neg' to make it easier to pull out only capacity
    % transients of interest without having to check the notebook.
    protName = 'cc_gapfree';
    % two alternatives for finding instances of the desired protocol
    % find(~cellfun('isempty',strfind(ephysData.(cellName).protocols,'ct_neg')));
    protLoc = find(strncmp(protName,ephysData.(cellName).protocols,6));
    
    if protLoc
        for i = 1:length(protLoc)
            gapfree(:,i) = ephysData.(cellName).data{1,protLoc(i)};
        end
        basalVoltage{iCell} = gapfree;
    end
end

clear allCells iCell cellName protName protLoc i gapfree;