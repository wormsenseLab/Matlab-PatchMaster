protList ={'_time'};
matchType = 'last';
strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 150];
extFilterFreq = 5;
includeFlag = 1;

antSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist,'included', 1, ...
     'stimFilterFrequencykHz', extFilterFreq);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%%
ExcludeSweeps(ephysData, protList, antSineCells, 'matchType', matchType, 'channel', 1);

clear protList matchType;

%%
sinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

% okay we're just going to plot the power spectrum in the morning and
% ignore the bode plot because I don't know what a "system" object really
% is in matlab or if my data can be one, or what a linear time-invariant
% system is or why a Bode plot is useful
%%
eachFreq = [0 10 30 100 200 500 1000];
sf = 10000; %Hz
allPSD = [];
allF = [];
allSize = [];
% ally = []; %for plotting fft magnitude and phase
% allf_fft = [];


for iRec = 1:size(sinePeaksNorm,1)
    theseRecs = sinePeaksNorm{iRec,2};
    theseSizes = sinePeaksNorm{iRec,3}(:,1:3);
    
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

%% Plot steady state sine comparison

eachFreq = [0 10 30 100 200 500 1000];
sf = 10000; %Hz
allSteady = [];
allSize = [];

for iRec = 1:size(sinePeaksNorm,1)
    theseSizes = sinePeaksNorm{iRec,3}(:,[1:3,5]);
    
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
