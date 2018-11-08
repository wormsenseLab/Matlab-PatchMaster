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
% [ephysData,tree] = ImportPatchData();
ephysData = ImportPatchData('incl',1);
% ephysData = ImportPatchData(ephysData, 'incl',1);


% Keep only data with given project prefixes/names.
projects = {'FAT';'SYM'};
% projects = {'SAR','FAT'};

%fix a couple messed up names
ephysData.FAT029 = ephysData.FAT029s;
ephysData.FAT017 = ephysData.FAT017e001;
ephysData.FAT164 = ephysData.FAT164001;
ephysData = rmfield(ephysData,{'FAT017e001';'FAT029s';'FAT164001'});

%get rid of the crap/unrelated recordings
ephysData = FilterProjectData(ephysData, projects);

clear projects;
%% Analyze capacity transient for C, Rs, and tau

ephysData = CtAnalysis(ephysData);
% ephysData = CtAnalysis(ephysData,newCells);

clear lastfit;

%% Import metadata with info about which IVq protocols to look at

ephysIVMetaData = ImportMetaData(); %FAT-IV Assignments

%% Import Recording Database
ephysMetaData = ImportMetaData();  %Recording Database

%% Make list of approved traces (by selecting traces to exclude)

% protList = {'DispRate'};
% protList = {'PrePulse'};
% protList = {'WC_Probe';'WC_ProbeSmall';'WC_ProbeLarge'};
protList ={'WC_Probe';'NoPre'};

matchType = 'first';
strainList = {'TU2769'};
internalList = {'IC2'};
% cellTypeList = {'ALMR'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 100];
resistCutoff = '<250';
extFilterFreq = [2.5 5];
includeFlag = 1;

anteriorStepCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
     'stimFilterFrequencykHz', extFilterFreq, 'included', 1);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%% Run sweep selection GUI

ExcludeSweeps(ephysData, protList, anteriorStepCells, 'matchType', matchType);

clear protList matchType;

%% Generic IdAnalysis run

% protList ={'DispRate'};
% sortSweeps = {'velocity','magnitude','magnitude','magnitude'};
% matchType = 'first';
% rateMRCs = IdAnalysis(ephysData,protList,rateCells,'num','matchType',matchType,...
%     'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1, ...
%     'saveSweeps',1);


% protList = {'Pair8'};
% matchType = 'last';
% sortSweeps = {'magnitude','magnitude','interval','magnitude'};
% testMRCs = IdAnalysis(ephysData,protList,matchType,'num','sortSweepsBy',sortSweeps);


% protList ={'PrePulse'};
% sortSweeps = {'position','position','position','position'};
% matchType = 'full';
% wtPreMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
%     'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1);

% protList ={'Noise_Trap'};
% sortSweeps = {'velocity','velocity','magnitude','magnitude'};
% matchType = 'first';
% rateMRCs = IdAnalysis(ephysData,protList,preCells,'num','matchType',matchType,...
%     'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1);


% protList ={'WC_Probe8'};
protList = {'WC_Probe','NoPre'};
sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};
matchType = 'first';
anteriorMRCs = IdAnalysis(ephysData,protList,anteriorStepCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1);
clear protList sortSweeps matchType

%% NonStat Noise Analysis


protList ={'NoisePre'};
matchType = 'first';
noisePre = NonStatNoiseAnalysis(ephysData,protList,preCells,'matchType',matchType);
clear protList matchType


%% Frequency Dependence Analysis

protList ={'_time'};
matchType = 'last';
sinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);
clear protList matchType

allSines = vertcat(sinePeaks{:,3});
sortedSines = sortrows(allSines, 3);
%% Print list of Rs

recs = fieldnames(ephysData);
resists = nan(length(recs),3);

for i=1:length(recs)
try fprintf('%s: %s\n', recs{i}, sprintf('%6g',round(ephysData.(recs{i}).Rs(ephysData.(recs{i}).protRs))))
catch
    continue
end
resists(i,1) = str2double(recs{i}(end-2:end));
resists(i,2) = ephysData.(recs{i}).C(ephysData.(recs{i}).protRs);
resists(i,3) = ephysData.(recs{i}).Rs(ephysData.(recs{i}).protRs);

end

clear i 

%% Print all Rs for making FAT_IV_Assignments spreadsheet

%prints list with recording name and three rows: C, Rs, and series number
%(use when deciding not only which recordings, but which protocols meet the
%series resistance requirement)

recs = fieldnames(ephysData);
% recs = newCells;

for i=1:length(recs)
try fprintf('%s: %s\n', recs{i}, sprintf('%6g',round(ephysData.(recs{i}).C * 1e12,2)))
catch
    continue
end
fprintf('%s %s\n', '       ', sprintf('%6g',round(ephysData.(recs{i}).Rs)))
fprintf('%s %s\n', '       ', sprintf('%6g',ephysData.(recs{i}).RsSeries));
end

%% Assign numbers of IVq series to look at

% Define anonymous function to convert cells from metadata to double arrays
% for protStart and protRs
CellToArray = @(x) reshape([x{:}],size(x,1),size(x,2), size(x,3));

% List which sets of ct_ivqs to use for on cell (col 2)/whole cell (col 3)
% calculations for the above selection of cells. (i.e., if the second  
% ct_ivq protocol run for that recording was the one you want to use for  
% whole-cell, put "2" in col 2 for that recording). Make sure it has three 
% sequential ivq pgfs.
allCells = ephysIVMetaData(2:end,1);

%TODO: don't assume headers: try/catch it
protStart = ephysIVMetaData(2:end,2:3)'; 
protStart = CellToArray(protStart);

protRs = ephysIVMetaData(2:end,4)';
protRs = CellToArray(protRs);


for i = 1:length(allCells)
    ephysData.(allCells{i}).protOC = protStart(1,i);
    ephysData.(allCells{i}).protWC = protStart(2,i);    
    ephysData.(allCells{i}).protRs = protRs(i);
end
clear i protStart protRs CellToArray
%% Process voltage steps

% allCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT027';'FAT028';'FAT029';'FAT030';'FAT031';'FAT032'};
% 
% allIVs = IVAnalysis(ephysData,allCells);

% wtCells = {'FAT020';'FAT021';'FAT022';'FAT025';'FAT031';'FAT033';'FAT034';'FAT035'};
% fatCells = {'FAT027';'FAT028';'FAT029'; 'FAT030'; 'FAT032';
%     'FAT036';'FAT037';'FAT038';'FAT039';'FAT040';'FAT041';'FAT042';
%     'FAT043';'FAT044'};
% 
% testingSplit = IVAnalysis(ephysData,wtCells,fatCells);
% wtIVs = testingSplit{1};
% fatIVs = testingSplit{2};
% 
% clear testingSplit wtCells fatCells

allIVs = IVAnalysis(ephysData,allCells);
allIVs = IVRsCorrection(ephysData,allIVs,allCells);

allI = nan(1050,12,length(allCells));
allV = nan(12,length(allCells));
whichCells = false(length(allCells),1);
for i = 1:length(allCells)
   try allI(:,:,i) = allIVs.(allCells{i}).capCorrIV;
   catch
       continue
   end
   try allV(:,i) = allIVs.(allCells{i}).actualV;
   catch 
       continue
   end
   whichCells(i) = true;
end
ivCells = allCells(whichCells);

genotype = cell(length(ivCells),2);
for i=1:length(ivCells)
genotype(i,1) = ivCells(i);
genotype(i,2) = ephysMetaData(strcmp(ephysMetaData(:,1),ivCells(i)),2);
end
wtCells = ivCells(strcmp(genotype(:,2),'TU2769'));
fatCells = ivCells(strcmp(genotype(:,2),'GN381'));

wtI = allI(:,:,ismember(ivCells,wtCells));
fatI = allI(:,:,ismember(ivCells,fatCells));
wtV = allV(:,ismember(ivCells,wtCells));
fatV = allV(:,ismember(ivCells,fatCells));

wtIsteady = squeeze(mean(wtI(350:550,:,:)));
fatIsteady = squeeze(mean(fatI(350:550,:,:)));

meanwtI = nanmean(wtIsteady,2);
meanwtV = nanmean(wtV,2);
steWTI = nanstd(wtIsteady,0,2)./sqrt(sum(~isnan(wtIsteady(1,:))));
steWTV = nanstd(wtV,0,2)./sqrt(sum(~isnan(wtV(1,:))));

meanfatI = nanmean(fatIsteady,2);
meanfatV = nanmean(fatV,2);
steFatI = nanstd(fatIsteady,0,2)./sqrt(size(fatIsteady,2));
steFatV = nanstd(fatV,0,2)./sqrt(size(fatV,2));

% wtIVs = rmfield(allIVs, allCells(~ismember(allCells,wtCells)));
% fatIVs = rmfield(allIVs, allCells(~ismember(allCells,fatCells)));
%% Plot mechanically evoked currents in response to single steps
% IdAnalysis

% wtCells = {'FAT034';'FAT035'};
% fatCells = {'FAT036';'FAT038';'FAT042';'FAT043';'FAT044'};
% % FAT041 has high leak

mechPeaksWT = IdAnalysis(ephysData,wtCells,1);
mechPeaksFat = IdAnalysis(ephysData,fatCells,1);

% mechPeaks = IdAnalysis(ephysData,allCells,0);
mechCellsWT = allCells(~cellfun('isempty',mechPeaksWT(:,1)));
mechPeaksWT = mechPeaksWT(~cellfun('isempty',mechPeaksWT(:,1)),:);
mechCellsFat = allCells(~cellfun('isempty',mechPeaksFat(:,1)));
mechPeaksFat = mechPeaksFat(~cellfun('isempty',mechPeaksFat(:,1)),:);


%% Look at interstimulus interval
% allCells = {'FAT059'; 'FAT061';'FAT062';'FAT063'};
%     
% ISIs = ISIAnalysis(ephysData,allCells);
% 
int1sCells = {'FAT059'; 'FAT061';'FAT062';'FAT063';'FAT064'};
int3sCells = {'FAT065';'FAT066';'FAT072';'FAT073';'FAT077'};

isi1s = ISIAnalysis(ephysData,int1sCells,'WC_Probe8');
isi3s = ISIAnalysis(ephysData,int3sCells,'WC_Probe8_3s');

for i = 1:5
    clear a
    a(1,:)=isi1s{3,i}(1):1:isi1s{3,i}(1)+length(isi1s{1,i})-1;
    if size(isi1s{3,i},2)>1
        a(2,:)=isi1s{3,i}(2):1:isi1s{3,i}(2)+length(isi1s{1,i})-1;
    end
    tVec1{i}=reshape(a',[],1);
    iVecOn1{i} = reshape(isi1s{4,i}',[],1);
    iVecOff1{i} = reshape(isi1s{5,i}',[],1);
end

for i = 1:5
    clear a
    a(1,:)=isi3s{3,i}(1):1:isi3s{3,i}(1)+length(isi3s{1,i})-1;
    if size(isi3s{3,i},2)>1
        a(2,:)=isi3s{3,i}(2):1:isi3s{3,i}(2)+length(isi3s{1,i})-1;
    end
    tVec3{i}=reshape(a',[],1);
    iVecOn3{i} = reshape(isi3s{4,i}',[],1);
    iVecOff3{i} = reshape(isi3s{5,i}',[],1);
end

figure(); 
h1s = subplot(1,2,1); hold on;
h3s = subplot(1,2,2); hold on;

for i=1:5
    plot(h1s,tVec1{i},iVecOn1{i},'b');
    plot(h1s,tVec1{i},iVecOff1{i},'r');
    plot(h3s,tVec3{i},iVecOn3{i},'b');
    plot(h3s,tVec3{i},iVecOff3{i},'r');    
end

%% Read in genotypes

allCells = fieldnames(ephysData);

genotype = cell(length(allCells),2);
for i=1:length(allCells)
genotype(i,1) = allCells(i);
try genotype(i,2) = ephysMetaData(strcmp(ephysMetaData(:,1),allCells(i)),2);
catch
    continue
end
end

clear i

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


%% Dt 
allCells = {'FAT032'};

OnDtTest = OnDtAnalysis(ephysData,allCells);

figure()
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.on1))
hold on;
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.on2),'d')
scatter(OnDtTest.FAT032.dts,mean(OnDtTest.FAT032.off),'s')


%% PrePulse/NoPrePulse
protList = {'PrePulse'};
matchType = 'first';
prepulseMRCs = IdAnalysis(ephysData,0, 'time',protList, matchType);

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