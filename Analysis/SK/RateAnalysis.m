% RateAnalysis.m
%
%

rateTracePicks = ImportMetaData();
rateTracePicks = metaDataConvert(rateTracePicks);

% stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
stimComI = ephysData.FAT075.data{2,11} ./ 0.408;
indentI = -ephysData.FAT075.data{3,11};
threshFraction = 0.05;

[~,rampEnd] = findSteps(5,stimComI,5,0.1,'endTime',300);
stepPDSize = findSteps(5,indentI,5,0.3,'endTime',300); %find this thresh based on noise, like in findMRCs
[~,rampEnd] = findSteps(5,indentI,5,stepPDSize(5)*.5,'endTime',300);

[~,rampEnd] = findSteps(5,stimComI,5,stepComSize(1)*threshFraction,'endTime',300);
[~, ~, rampStart] = findSteps(5,stimComI,5,stepComSize(1)*threshFraction);



threshTime = 30; % use first n ms of trace for setting noise threshold
sf=5;
baseTime = 30;

for iSeries = 1:length(rateTracePicks)
    cellName = rateTracePicks{iSeries,1};
    thisSeries = rateTracePicks{iSeries,1};
    
    stimComI = ephysData.(cellName).data{2,thisSeries} ./ 0.408;
    indentI = -ephysData.(cellName).data{3,thisSeries};
    probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
    leakSubtract = ...
        SubtractLeak(probeI, sf, 'BaseLength', baseTime);
    
    
    nSweeps = length(rateTracePicks{iSeries,3});
    
    rampThresh = 1.5*thselect(stimComI(1:threshTime*sf),'rigrsure');
    %stimWindow MUST BE VECTOR
    stimWindow(1) = 750;
    [stimSize,stimWindow(2)] = findSteps(nSweeps,stimComI,sf,rampThresh,'endTime',300);
    
        
        
end
