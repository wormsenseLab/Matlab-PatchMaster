strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 200];
extFilterFreq = [1 5];
includeFlag = 1;
cultT = 15;

wtSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist,'included', 1, ...
     'stimFilterFrequencykHz', extFilterFreq, 'TCultC', cultT);
 
 strainList = {'GN381'};
 
 fatSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist,'included', 1, ...
     'stimFilterFrequencykHz', extFilterFreq, 'TCultC', cultT);


clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep includeFlag extFilterFreq;

%NEXT: check to see if there's a difference in the power spectrum (and
%amplitude attenuation where applicable) for 2.5 vs. 5kHz photodiode
%signals. If they look the same, combine the data.

%%
protList ={'_time'};
matchType = 'last';

ExcludeSweeps(ephysData, protList, wtSineCells, 'matchType', matchType, 'channel', 1);
ExcludeSweeps(ephysData, protList, fatSineCells, 'matchType', matchType, 'channel', 1);

% protList ={'_num'};
% ExcludeSweeps(ephysData, protList, fatSineCells, 'matchType', matchType, 'channel', 1);

clear matchType;

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

protList ={'_time'};
matchType = 'last';

wtSinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

wtSinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0);

fatSinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

fatSinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0);


protList ={'_num'};
matchType = 'last';

fatSinePeaksNum = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'multiRep',1,'steadyTime',50, 'sortStimBy','endTime');

fatSinePeaksNumNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'multiRep',1,'steadyTime',50, 'sortStimBy','endTime','norm',1);
% sinePeaksPD = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0, 'channel',3);
% sinePeaksStim = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0, 'channel',2);

% okay we're just going to plot the power spectrum in the morning and
% ignore the bode plot because I don't know what a "system" object really
% is in matlab or if my data can be one, or what a linear time-invariant
% system is or why a Bode plot is useful
%% Calculate PSD
figure();
sf = 10000; %Hz
allPSD = [];
allF = [];
allFreq = [];

whichPeaks = wtSinePeaks; %also change meansByFreqI at bottom
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
        [pxx,f] = pwelch(thisRec,500,[],[],sf);
%         y = fft(thisRec,length(thisRec));
%         f_fft = ((0:1/length(thisRec):1-1/(length(thisRec)))*sf).';
        if ~isempty(allF) && length(allF)>length(f)
            f = [f; nan(length(allF)-length(f),size(f,2))];
            pxx = [pxx; nan(length(allPSD)-length(pxx),size(pxx,2))];            
        end
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

eachFreq = unique(sortedFreq);
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
legend(num2str(eachFreq));
box off;
chH = get(gca,'children');
set(gca,'children',flipud(chH));
plotfixer;

wtMeanPSDByFreq = meanPSDByFreq;

%% Plot stacked power spectra
figure(); clear axh;
yyh = cell(0);

for i = 1:size(meanPSDByFreqI,2)
    axh(i) = subtightplot(size(meanPSDByFreqI,2),1,i,0.02,0.05,0.1);
    yyh{i} = plotyy(meanF,meanPSDByFreqI(:,i),meanF,meanPSDByFreqPD(:,i));
    set(yyh{i}, 'YScale', 'log', 'XScale', 'log','box','off');
    set(yyh{i}(1),'YLim',[1e-30 1e-20])
    set(yyh{i}(2),'YLim',[1e-10 1e0])
end

for i = 1:size(meanPSDByFreqI,2)-1
    set(yyh{i},'XTickLabel',[]);
end
    

%% Calculate steady state mean and rms by frequency

whichPeaks = wtSinePeaks;
normFlag = 1; % 1 to normalize steady-state and RMS to value for square pulse
extFilter = [2.5 5];

sf = 10000; %Hz
allSteady = [];
allRMS = [];
allFreq = [];
allExtFilt = [];
allSquare = [];

for iRec = 1:size(whichPeaks,1)
    theseFreqs = whichPeaks{iRec,3}(:,[1:4,6,7,8]);
    if normFlag
        try theseFreqs(:,5:7) = bsxfun(@rdivide,theseFreqs(:,5:7),theseFreqs(:,7));
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
eachFreq = unique(sortedFreq);


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

wtFreq = eachFreq;
wtSteady = meanSteadyByFreq;
wtRMS = meanRMSByFreq;

%% Plot steady states

figure(); hold on;
errorbar(wtFreq,wtSteady(:,1),wtSteady(:,3));
errorbar(fatFreq,fatSteady(:,1),fatSteady(:,3));

figure(); hold on;
errorbar(wtFreq,wtRMS(:,1),wtRMS(:,3));
errorbar(fatFreq,fatRMS(:,1),fatRMS(:,3));

%% Export everything for Igor plotting and source data

