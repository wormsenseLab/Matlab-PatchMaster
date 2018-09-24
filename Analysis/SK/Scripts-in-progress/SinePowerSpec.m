protList ={'_time'};
matchType = 'last';
strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 150];
extFilterFreq = [2.5 5];
includeFlag = 1;

antSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist,'included', 1, ...
     'stimFilterFrequencykHz', extFilterFreq);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

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
selectedSweeps(:,3) = cellfun(@(x) ['' x],selectedSweeps(:,3));

% save the new selectedSweeps list to file
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save sweep list to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName, selectedSweeps);

clear fName filename pathname a b selectedSweeps_I selectedSweeps_PD selectedSweeps selectIdx


%%
sinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

sinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0);

sinePeaksPD = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 0, 'channel',3);

% okay we're just going to plot the power spectrum in the morning and
% ignore the bode plot because I don't know what a "system" object really
% is in matlab or if my data can be one, or what a linear time-invariant
% system is or why a Bode plot is useful
%%
figure();
eachFreq = [0 10 30 100 200 500 1000];
sf = 10000; %Hz
allPSD = [];
allF = [];
allSize = [];

whichPeaks = sinePeaksPD;
% ally = []; %for plotting fft magnitude and phase
% allf_fft = [];


for iRec = 1:size(whichPeaks,1)
    theseRecs = whichPeaks{iRec,2};
    theseSizes = whichPeaks{iRec,3}(:,1:3);
    
    for iFreq = 1:length(theseRecs)
        try thisRec = theseRecs{iFreq}(theseSizes(iFreq,1):theseSizes(iFreq,2),:);
        %use pwelch instead of periodogram to allow for use of Hamming
        %window (reduces variance but decreases peak resolution/increases peak width)
        %alternatively, periodogram with shorter window (or full-length
        %hamming window to reduce edge effects?)
        [pxx,f] = pwelch(thisRec,5000,[],[],sf);
%         y = fft(thisRec,length(thisRec));
%         f_fft = ((0:1/length(thisRec):1-1/(length(thisRec)))*sf).';
        allF = [allF mean(f,2)];
        allPSD = [allPSD mean(pxx,2)];
        allSize = [allSize theseSizes(iFreq,3)];
%         ally = [ally mean(y,2)];
%         allf_fft = [allf_fft f_fft];
        catch
        end
    end
    
end

[sortedSize, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
    sortRowsTol(allSize', 1, 1);

sortedPSD = allPSD(:,sortIdx);
% sortedy = ally(:,sortIdx);
meanF = allF(:,1);
% meanf_fft = allf_fft(:,1);
groupIdx = cell(0);
meansByFreq = [];
% yMean = [];

for iProfile = 1:length(eachFreq)
    groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
    thesePSD = sortedPSD(:,groupIdx{iProfile});
    meansByFreq (:,iProfile) = mean(thesePSD,2);
%     yMean(:,iProfile) = mean(sortedy(:,groupIdx{iProfile}),2);
end
% 
% magnitudeY = abs(yMean(:,7));        % Magnitude of the FFT
% phaseY = unwrap(angle(yMean(:,7)));  % Phase of the FFT
% helperFrequencyAnalysisPlot1(meanf_fft,magnitudeY,phaseY,size(y,1))

plot(meanF,meansByFreq);
set(gca, 'YScale', 'log', 'XScale', 'log');

xlabel('frequency (Hz)');
ylabel('A^2');
legend(num2str(eachFreq'));
box off;
chH = get(gca,'children');
set(gca,'children',flipud(chH));
plotfixer;

% meansByFreqPD = meansByFreq;
% meansByFreqI = meansByFreq;

%% Plot stacked power spectra
figure();
yyh = cell(0);

for i = 1:size(meansByFreq,2)
    axh(i) = subtightplot(size(meansByFreq,2),1,i,0.02,0.05,0.1);
    yyh{i} = plotyy(meanF,meansByFreqI(:,i),meanF,meansByFreqPD(:,i));
    set(yyh{i}, 'YScale', 'log', 'XScale', 'log','box','off');
    set(yyh{i}(1),'YLim',[1e-30 1e-20])
    set(yyh{i}(2),'YLim',[1e-10 1e0])
end

for i = 1:size(meansByFreq,2)-1
    set(yyh{i},'XTickLabel',[]);
end
    

%% Plot steady state sine comparison

whichPeaks = sinePeaksNorm;

eachFreq = [0 10 30 100 200 500 1000];
sf = 10000; %Hz
allSteady = [];
allSize = [];

for iRec = 1:size(whichPeaks,1)
    theseSizes = whichPeaks{iRec,3}(:,[1:3,6]);
    
    for iFreq = 1:size(theseSizes,1)
        allSteady = [allSteady; theseSizes(iFreq,4)];
        allSize = [allSize; theseSizes(iFreq,3)];
        
        
    end
end


[sortedSize, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
    sortRowsTol(allSize, 1, 1);

sortedSteady = allSteady(sortIdx);
groupIdx = cell(0);
meansByFreq = [];
sdByFreq = [];
nByFreq = [];
semByFreq = [];

for iProfile = 1:length(eachFreq)
    groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
    theseSteady = sortedSteady(groupIdx{iProfile});
    meansByFreq (iProfile,1) = mean(theseSteady);
    sdByFreq (iProfile,1) = std(theseSteady);
    nByFreq (iProfile,1) = length(theseSteady);
    semByFreq (iProfile,1) = std(theseSteady)/sqrt(length(theseSteady));
end
