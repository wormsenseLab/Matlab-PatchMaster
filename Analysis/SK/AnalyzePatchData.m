% AnalyzePatchData.m
% 
% TODO: Figure out import for RecordingDatabase.xlsx, add columns for OC,WC
% series numbers of interest and pull those.
% TODO: Pass in series numbers from spreadsheet for analysis, rather than
% analyzing all series with that protocol name (but check protocol name of
% those series and pass error if not).

%% Import Data
%this is a test

% Don't forget to run sigTOOL first!
ephysData = ImportPatchData();

% Keep only data with given project prefixes/names.
projects = {'FAT'};

ephysData = FilterProjectData(ephysData, projects);

ephysMetaData = ImportMetaData();

clear projects
%% Analyze capacity transient for C, Rs, and tau

ephysData = CtAnalysis(ephysData);

%% Temporary: Assign numbers of IVq series to look at

% TODO: Read these in from a separate file.

allCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT027';'FAT028';
    'FAT029';'FAT030';'FAT031'; 'FAT032'; 'FAT033'; 'FAT034';'FAT035';
    'FAT036';'FAT037';'FAT038';'FAT039';'FAT040';'FAT041';'FAT042';
    'FAT043';'FAT044'};


% List which sets of ct_ivqs to use for on cell (row 1)/whole cell (row 2)
% calculations for the above selection of cells. (i.e., if the second  
% ct_ivq protocol run for that recording was the one you want to use for  
% whole-cell, put "2" in row 1). Make sure it has three sequential ivq
% pgfs.
protStart = [1 1 1 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; ...
             4 4 4 7 4 7 6 4 4 4 4 4 4 4 4 4 4 4 4 4 4 1];
         
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

wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031';'FAT033';'FAT034';'FAT035'};
fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032';
    'FAT036';'FAT037';'FAT038';'FAT039';'FAT040';'FAT041';'FAT042';
    'FAT043';'FAT044'};

testingSplit = IVAnalysis(ephysData,wtCells,fatCells);
wtIVs = testingSplit{1};
fatIVs = testingSplit{2};

clear testingSplit wtCells fatCells

%% Correct voltage steps based on series resistance and plot as scatter

wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031'; 'FAT033';'FAT034';'FAT035'};
fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032';
    'FAT036';'FAT037';'FAT038';'FAT039';'FAT040';'FAT041';'FAT042';
    'FAT043'};
% FAT044 has no OC

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

for i = 1:length(fatCells)
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
clear i cellName
%% Plot mechanically evoked currents in response to single steps
% IdAnalysis

wtCells = {'FAT034';'FAT035'};
fatCells = {'FAT036';'FAT038';'FAT042';'FAT043';'FAT044'};
% FAT041 has high leak

mechPeaksWT = IdAnalysis(ephysData,wtCells);
mechPeaksFat = IdAnalysis(ephysData,fatCells);


%% Look at interstimulus interval
allCells = {'FAT059'; 'FAT061';'FAT062';'FAT063'};
    
ISIs = ISIAnalysis(ephysData,allCells);


%% Plot single MRC sets
% Draw stim protocol for MRCs
dt = 0.2; % ms, based on sampling frequency (5kHz in current WC_Probe)
tVec = 0:dt:500-dt;
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
toPlot = ephysData.FAT029.data{1,11}(:,6);
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
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.on1))
hold on;
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.on2),'d')
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.off),'s')


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


%% Check linearity of photodiode response

% % Import file '150702_photodiodetest1.dat', then:
% testMed = ephysData.E1.data{3,43};
% testSml = ephysData.E1.data{3,45};
% testLrg = ephysData.E1.data{3,44};

testMedBase = mean(testMed(1:500,:));
testMedStep = mean(testMed(1000:2000,:));
testMedSub = testMedStep-testMedBase;

testSmlBase = mean(testSml(1:500,:));
testSmlStep = mean(testSml(1000:2000,:));
testSmlSub = testSmlStep-testSmlBase;

testLrgBase = mean(testLrg(1:500,:));
testLrgStep = mean(testLrg(1000:2000,:));
testLrgSub = testLrgStep-testLrgBase;

controlStim = [1 3 5 7 9 11 0.5 1.5 8 10];
pdResp = [testMedSub testSmlSub testLrgSub];
toExcel = [controlStim; pdResp]';

scatter(controlStim, pdResp);