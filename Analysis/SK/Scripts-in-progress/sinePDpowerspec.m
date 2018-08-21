protList ={'_time'};
matchType = 'last';
strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 150];
extFilterFreq = 5;

antSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, ...
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

for iRec = 1:size(sinePeaksNorm,1)
    theseRecs = sinePeaksNorm{iRec,2};
    theseSizes = sinePeaksNorm{iRec,3}(:,1:3);
    
    for iFreq = 1:length(theseRecs)
        try thisRec = theseRecs{iFreq}(theseSizes(iFreq,1):theseSizes(iFreq,2),:);
        [pxx,f] = periodogram(thisRec,[],[],sf);
        allF = [allF mean(f,2)];
        allPSD = [allPSD mean(pxx,2)];
        allSize = [allSize theseSizes(iFreq,3)];
        catch
        end
    end
    
end

[sortedSize, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
    sortRowsTol(allSize', 1, 1);

sortedPSD = allPSD(:,sortIdx);
meanF = allF(:,1);
groupIdx = cell(0);
meansByFreq = [];

for iProfile = 1:length(eachFreq)
    groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
    thesePSD = sortedPSD(:,groupIdx{iProfile});
    meansByFreq (:,iProfile) = mean(thesePSD,2);
end

plot(meanF,10*log10(meansByFreq));
xlabel('frequency (Hz)');
ylabel('dB');
legend(num2str(eachFreq'));
box off;
% plotfixer;

chH = get(gca,'children');
set(gca,'children',flipud(chH));