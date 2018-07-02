protList ={'_time'};
matchType = 'last';
strainList = {'TU2769'};
internalList = {'IC2'};
stimPosition = {'anterior'};
wormPrep = {'dissected'};
cellDist = [40 150];
extFilterFreq = 2.5;

antSineCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, ...
     'stimFilterFrequencykHz', extFilterFreq);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%%
ExcludeSweeps(ephysData, protList, antSineCells, 'matchType', matchType, 'channel', 3);

clear protList matchType;

%%
sinePeaksNorm = FrequencyAnalysis(ephysData, ephysMetaData, protList, 'matchType', matchType, 'norm', 1);

