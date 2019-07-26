% SinePowerSpec.m
% 
% Analysis of responses to sinusoidal stimuli. Outputs plots of power
% spectral density for each frequency for both the current response and the
% stimulus photodiode signal. Calculates mean steady-state current and rms 
% at the end of the sine and plots against simulated data.

%%Import data and metadata
ephysData = ImportPatchData('incl',1);
projects = {'FAT';'SYM'};
ephysData = FilterProjectData(ephysData, projects);
%fix a couple messed up names
ephysData.FAT029 = ephysData.FAT029s;
ephysData.FAT017 = ephysData.FAT017e001;
ephysData.FAT164 = ephysData.FAT164001;
ephysData = rmfield(ephysData,{'FAT017e001';'FAT029s';'FAT164001'});

clear projects;

ephysMetaData = ImportMetaData();  %Recording Database.xlsx
attenuationData = ImportMetaData(); %AttenuationCalcs_IC2_190620.xlsx

%% Select by recording conditions
protList ={'_time'};
matchType = 'last';
strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 200];
extFilterFreq = [2.5 5];
resistCutoff = '<250'; % Rs < 250 MOhm
includeFlag = 1;

antSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist,'included', includeFlag, ...
     'stimFilterFrequencykHz', extFilterFreq,'RsM',resistCutoff);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep includeFlag extFilterFreq;

%NEXT: check to see if there's a difference in the power spectrum (and
%amplitude attenuation where applicable) for 2.5 vs. 5kHz photodiode
%signals. If they look the same, combine the data.

%%

ExcludeSweeps(ephysData, protList, antSineCells, 'matchType', matchType, 'channel', 1);

clear protList matchType;

%% Check if the PD signal is useful/shows abnormalities
% Make list of only sweeps with stable current AND good PD signal (lower
% noise, high enough voltage for SNR, no weird movements during sine)

% selectedSweeps_I = ExcludeSweeps(ephysData, protList, antSineCells, 'matchType', matchType, 'channel', 1);
% selectedSweeps_PD = ExcludeSweeps(ephysData, protList, antSineCells, 'matchType', matchType, 'channel', 3);

selectedSweeps_I = ImportMetaData();
selectedSweeps_PD = ImportMetaData();


% This method of comparing the two lists only works if there's only one
% series/sweep per row. It simply concatenates cell ID, series, and sweep
% number into a char row. For more sweeps per series, you will need to loop
% (find matching cell/series combo, then check which sweeps are matching).
a = selectedSweeps_I;
b = selectedSweeps_PD;
a(:,2) = cellfun(@num2str, a(:,2),'un',0);
b(:,2) = cellfun(@num2str, b(:,2),'un',0);

a = cellfun(@(x,y,z) [x y z], a(:,1),a(:,2),a(:,3),'un',0);
b = cellfun(@(x,y,z) [x y z], b(:,1),b(:,2),b(:,3),'un',0);
a = vertcat(a{:});
b = vertcat(b{:});

selectIdx = ismember(a,b,'rows');

selectedSweeps = selectedSweeps_I(selectIdx,:);
selectedSweeps(:,3) = cellfun(@(x) ['' x],selectedSweeps(:,3),'un',0);

% save the new selectedSweeps list to file
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save sweep list to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName, selectedSweeps);

clear fName filename pathname a b selectedSweeps_I selectedSweeps_PD selectedSweeps selectIdx


%%
% sinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

sinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList,antSineCells, 'matchType', matchType, 'norm', 0);

sinePeaksPD = FrequencyAnalysis(ephysData, ephysMetaData, protList,antSineCells, 'matchType', matchType, 'norm', 0, 'channel',3);
sinePeaksStim = FrequencyAnalysis(ephysData, ephysMetaData, protList,antSineCells, 'matchType', matchType, 'norm', 0, 'channel',2);

% sine_allExt_ant_PDfiltered(180923).xls for all

%% Calculate PSD
figure();
eachFreq = [0 10 30 100 200 500 1000];
sf = 10000; %Hz
allPSD = [];
allF = [];
allFreq = [];

whichPeaks = sinePeaks; %also change meansByFreqI at bottom
% ally = []; %for plotting fft magnitude and phase
% allf_fft = [];


for iRec = 1:size(whichPeaks,1)
    theseRecs = whichPeaks{iRec,2};
    theseFreqs = whichPeaks{iRec,3}(:,1:3);
    
    for iFreq = 1:length(theseRecs)
        try thisRec = theseRecs{iFreq}(theseFreqs(iFreq,1):theseFreqs(iFreq,2),:);
        %use pwelch instead of periodogram to allow for use of Hamming
        %window (reduces variance but decreases peak resolution/increases peak width)
        %alternatively, periodogram with shorter window (or full-length
        %hamming window to reduce edge effects?)
        [pxx,f] = pwelch(thisRec,5000,[],[],sf);
%         y = fft(thisRec,length(thisRec));
%         f_fft = ((0:1/length(thisRec):1-1/(length(thisRec)))*sf).';
        allF = [allF mean(f,2)];
        allPSD = [allPSD mean(pxx,2)];
        allFreq = [allFreq theseFreqs(iFreq,3)];
%         ally = [ally mean(y,2)];
%         allf_fft = [allf_fft f_fft];
        catch
        end
    end
    
end

[sortedFreq, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
    sortRowsTol(allFreq', 1, 1);

sortedPSD = allPSD(:,sortIdx);
% sortedy = ally(:,sortIdx);
meanF = allF(:,1);
% meanf_fft = allf_fft(:,1);
groupIdx = cell(0);
meanPSDByFreq = [];
% yMean = [];

for iProfile = 1:length(eachFreq)
    groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
    thesePSD = sortedPSD(:,groupIdx{iProfile});
    meanPSDByFreq (:,iProfile) = mean(thesePSD,2);
%     yMean(:,iProfile) = mean(sortedy(:,groupIdx{iProfile}),2);
end
% 
% magnitudeY = abs(yMean(:,7));        % Magnitude of the FFT
% phaseY = unwrap(angle(yMean(:,7)));  % Phase of the FFT
% helperFrequencyAnalysisPlot1(meanf_fft,magnitudeY,phaseY,size(y,1))

plot(meanF,meanPSDByFreq);
set(gca, 'YScale', 'log', 'XScale', 'log');

xlabel('frequency (Hz)');
ylabel('A^2');
legend(num2str(eachFreq'));
box off;
chH = get(gca,'children');
set(gca,'children',flipud(chH));
plotfixer;

% meanPSDByFreqPD = meanPSDByFreq;
meanPSDByFreqI = meanPSDByFreq;

%% Read in simulated power spectra and steady-state/rms
freqSim_path = 'C:\Users\Sammy\Dropbox\Goodman Lab\Posters Papers Proposals\2018-Katta-SpatiotemporalDynamics\source-data\temporal-frequency\';

sim_path = fullfile(freqSim_path,'figure_4c_stimulus.txt');
sim_stim = dlmread(sim_path,'\t',1,0);
sim_f = sim_stim(:,1);
sim_stim = sim_stim(:,2:end);

sim_path = fullfile(freqSim_path,'figure_4c_response.txt');
sim_resp = dlmread(sim_path,'\t',1,1);

sim_path = fullfile(freqSim_path, 'figure_4d.txt');
sim_steady_t2p9 = dlmread(sim_path,'\t',1,0);
sim_sum_f = sim_steady_t2p9(:,1);
sim_steady_t2p9 = sim_steady_t2p9(:,2:end);

sim_path = fullfile(freqSim_path, 'figure_4e.txt');
sim_rms_t2p9 = dlmread(sim_path,'\t',1,1);

sim_path = fullfile(freqSim_path, 'figure_7d.txt');
sim_steady_t0p6 = dlmread(sim_path,'\t',1,0);
sim_sum_f_t0p6 = sim_steady_t0p6(:,1);
sim_steady_t0p6 = sim_steady_t0p6(:,2:end);

sim_path = fullfile(freqSim_path, 'figure_7e.txt');
sim_rms_t0p6 = dlmread(sim_path,'\t',1,1);

%% Plot stacked power spectra for sim and experimental

scrsz = get(groot,'ScreenSize');

% Plot experimental after removing 1kHz 
meanPSDByFreqPD = meanPSDByFreqPD(:,~ismember(eachFreq,1000));
meanPSDByFreqI = meanPSDByFreqI(:,~ismember(eachFreq,1000));
eachFreq = eachFreq(~ismember(eachFreq,1000));

figure('Position',[1 1 800 700]); clear axh;
yyh = cell(0);

for i = 1:size(meanPSDByFreqI,2)
    axh(i) = subtightplot(size(meanPSDByFreqI,2),1,i,0.02,0.05,0.1);
    yyh{i} = plotyy(meanF,meanPSDByFreqI(:,i),meanF,meanPSDByFreqPD(:,i),@loglog);
    set(yyh{i}(1),'YLim',[1e-30 1e-20])
    set(yyh{i}(2),'YLim',[1e-10 1e0])
    set(yyh{i}, 'YScale', 'log', 'XScale', 'log','box','off');
end

for i = 1:size(meanPSDByFreqI,2)-1
    set(yyh{i},'XTickLabel',[]);
end

% Plot simulated data
figure('Position',[1 1 800 700]); clear axh;
yyh = cell(0);

for i = 1:size(sim_resp,2)
    axh(i) = subtightplot(size(meanPSDByFreqI,2),1,i,0.02,0.05,0.1);
    yyh{i} = plotyy(sim_f,sim_resp(:,i),sim_f,sim_stim(:,i),@semilogx);
    set(yyh{i},'XLim',[1 1e5],'box','off');
    set(yyh{i}(1),'YLim',[-150 50])
    set(yyh{i}(2),'YLim',[-150 50])
end

for i = 1:size(sim_resp,2)-1
    set(yyh{i},'XTickLabel',[]);
end
    
plotfixer;
%% Calculate steady state mean and rms by frequency

whichPeaks = sinePeaks;
normFlag = 0; % 1 to normalize steady-state and RMS to value for square pulse
extFilter = [2.5 5];

eachFreq = [0 10 30 100 200 500];
sf = 10000; %Hz
allSteady = [];
allRMS = [];
allFreq = [];
allExtFilt = [];
allSquare = [];

for iRec = 1:size(whichPeaks,1)
    theseFreqs = whichPeaks{iRec,3}(:,[1:4,6,7,8]);
    if normFlag
        try theseFreqs(:,[5,6,7]) = bsxfun(@rdivide,theseFreqs(:,[5:7]),theseFreqs(theseFreqs(:,3)==10,[5:7]));
        catch
            continue
        end
    end
    for iFreq = 1:size(theseFreqs,1)
        allSteady = [allSteady; theseFreqs(iFreq,5)];
        allRMS = [allRMS; theseFreqs(iFreq,6)];
        allFreq = [allFreq; theseFreqs(iFreq,3)];
        allExtFilt = [allExtFilt; theseFreqs(iFreq,4)];
        allSquare = [allSquare; theseFreqs(iFreq,7)];
    end
end


[sortedFreq, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
    sortRowsTol(allFreq, 1, 1);

sortedSteady = allSteady(sortIdx);
sortedRMS = allRMS(sortIdx);
sortedSquare = allSquare(sortIdx);
sortedExtFilt = allExtFilt(sortIdx);
groupIdx = cell(0);
meanSteadyByFreq = [];
nByFreq = [];
meanRMSByFreq = [];


for iProfile = 1:length(eachFreq)
    groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
    whichFilt = ismember(sortedExtFilt(groupIdx{iProfile}),extFilter);
    if normFlag
        theseSteady = sortedSteady(groupIdx{iProfile}(whichFilt))./abs(sortedSquare(groupIdx{iProfile}(whichFilt)));
        theseRMS = sortedRMS(groupIdx{iProfile}(whichFilt))./abs(sortedSquare(groupIdx{iProfile}(whichFilt)));
    else
        theseSteady = sortedSteady(groupIdx{iProfile}(whichFilt));
        theseRMS = sortedRMS(groupIdx{iProfile}(whichFilt));
    end
    
    meanSteadyByFreq (iProfile,1) = mean(theseSteady);
    meanSteadyByFreq (iProfile,2) = std(theseSteady);
    meanSteadyByFreq (iProfile,3) = std(theseSteady)/sqrt(length(theseSteady));

    meanRMSByFreq (iProfile,1) = mean(theseRMS);
    meanRMSByFreq (iProfile,2) = std(theseRMS);
    meanRMSByFreq (iProfile,3) = std(theseRMS)/sqrt(length(theseRMS));

    nByFreq (iProfile,1) = length(theseSteady);
end


%% Plot experimental and simulated steady-state and RMS


figure(); clear axh;
yyh = cell(0);

yyh = plotyy(eachFreq,meanSteadyByFreq(:,1),sim_sum_f,sim_steady_t2p9,@semilogx);
set(yyh(1),'YLim',[0 60],'YTick',(0:20:60));
hold(yyh(1),'on');
errorbar(yyh(1),eachFreq,meanSteadyByFreq(:,1),meanSteadyByFreq(:,3),'bo');
hold(yyh(2),'on');
plot(yyh(2),sim_sum_f_t0p6,sim_steady_t0p6);
set(yyh,'box','off');
xlabel('Frequency (Hz)')
ylabel(yyh(1),'Steady-state current (pA)');
ylabel(yyh(2),'Normalized steady-state current');

figure(); clear axh;
yyh = cell(0);

yyh = plotyy(eachFreq,meanRMSByFreq(:,1),sim_sum_f,sim_rms_t2p9,@semilogx);
set(yyh(1),'YLim',[0 8],'YTick',(0:2:8));
hold(yyh(1),'on');
errorbar(yyh(1),eachFreq,meanRMSByFreq(:,1),meanRMSByFreq(:,3),'bo');
hold(yyh(2),'on');
plot(yyh(2),sim_sum_f_t0p6,sim_rms_t0p6);
set(yyh, 'box','off');
xlabel('Frequency (Hz)')
ylabel(yyh(1),'Steady-state RMS (pA)');
ylabel(yyh(2),'Normalized steady-state RMS');

%% Plot representative traces (FAT170) at 10, 100, 500 Hz
% 2,4,6 = 10,100,500Hz
% originally FAT218 bc it had 1kHz stim, but that only has 2-3 reps and 
% lower signal-to-noise

% check mean nReps for each recording
sineReps = cell2mat(cellfun(@(x) mean(x(:,5)),sinePeaks(:,3),'un',0));
% FAT170 has ~5-6 reps


% If you change filtering, make sure this still points to the right trace
respTraces = sinePeaks{3,2}([2 4 6]);
respTraces = cellfun(@(x) mean(x,2),respTraces,'un',0);
respTraces = [respTraces{:}];
respMean = mean(respTraces(16499:17499,:),1);

stimTraces = sinePeaksPD{3,2}([2,4,6]);
stimTraces = cellfun(@(x) mean(x,2),stimTraces,'un',0);
stimTraces = [stimTraces{:}];
tVec = ((1:length(stimTraces))/sf)';

% Plot calibrated photodiode trace as alternative
a = ephysData.FAT170.data{3,52}; % probe steps
a = a-mean(a(1:1499));
b = ephysData.FAT170.data{3,53}-ephysData.FAT170.data{3,54}; %pd steps minus worm only pd steps
b = b-mean(b(1:1499));
c = [0 2.5 5 7.5 10 12.5];

% calibCells = {'FAT170'};
% calibData = cell(0);
% calibData{1} = b;
% 
% clear handles;
% plotData = calibData{:,1};
% handles = selectCalibSteps(plotData,calibCells{1});
% stepIdx = handles.cursorPoints(:,1);
% stepIdx = reshape(stepIdx,2,[])';
% for j = 1:size(stepIdx,1)
%     stepValues(1,j) = mean(plotData(stepIdx(j,1):stepIdx(j,2)))';
% end
% close;      
% d = stepValues(1:6);
d = [0.000719859832312041,-0.196562355775423,-0.460558335247858,-0.755470591227155,-1.11678911522346,-1.51308976348689];
e = stimTraces;
pdTraces = interp1(-d,c,-e,'pchip');

yyh = plotyy(tVec,respTraces,tVec,pdTraces);
hline(respMean(1),'b'); hline(respMean(2),'r'); hline(respMean(3),'m');

%% Plot zoomed traces for 5 cycles
% Calculate time
% 10Hz = 500ms
% 100Hz = 50ms
% 500Hz = 10ms
cycleLims = [1.25 1.7500; 1.7000 1.7500; 1.7400 1.7500];

yyh = cell(0);
for i = 1:3
    figure();
    yyh{i} = plotyy(tVec,respTraces(:,i),tVec,pdTraces(:,i));
    hline(respMean(i),'k:');
    set(yyh{i},'XLim',cycleLims(i,:));
    set(yyh{i}(2),'YLim',[0 10]);
    set(yyh{i}(1),'Ylim',[respMean(i)-15e-12 respMean(i)+15e-12]);
end

%% Separate and write representative traces for channel-sim-distance figure
clear reps repStim
% FAT170 mean traces at 10, 100, 500 Hz 
% stimEnd = [16013; 17455; 17494];
stimEnd = [16505; 17404; 17484];

sf = 10; %kHz
boxTime = 500*sf; %ms
onBox = [stimEnd-boxTime stimEnd];
repTrace = cell2mat(cellfun(@(x) mean(x,2),sinePeaks{3,2}([2 4 6]),'un',0)');
stimTrace = cell2mat(cellfun(@(x) mean(x,2)/0.408,sinePeaksStim{3,2}([2 4 6]),'un',0)');
for i = 1:3
    reps(:,i) = repTrace(onBox(i,1):onBox(i,2),i);
    repStim(:,i) = stimTrace(onBox(i,1):onBox(i,2),i);
end
tVec = ((1:length(reps))-boxTime/2-1)'./(sf*1000);
plot(repStim)

%% Use steps to plot resonance in stim
a = cellfun(@(x) x{:}, sinePeaksPD(:,2),'un',0);
b = [a{:}];
c = mean(b,2);
sf = 10000;
onStim = c(1505:4500);
offStim = c(4505:7500);

plot(onStim);
% [pxx,f] = pwelch(onStim,2000,[],[],sf);
% plot(f,pxx)
% set(gca,'XScale','log','YScale','log');