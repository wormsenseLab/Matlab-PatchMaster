[~,gidx] = intersect(genotype(:,1),testMRCs(:,1));
rateGenotypes = genotype(gidx,:);
wtTauMRCs = testMRCs(strcmp('TU2769',rateGenotypes(:,2)),:);
for i=1:size(wtTauMRCs,1)
wtTauMRCs{i,5}=repmat(wtTauMRCs{i,1},size(wtTauMRCs{i,2},1),1);
end
for i=1:size(wtTauMRCs,1)
wtTauMRCs{i,2}(:,end+1:9000)=NaN;
end
wtTauNames=vertcat(wtTauMRCs{:,5});
wtTauStats = vertcat(wtTauMRCs{:,3});
wtTauMRCTrace = vertcat(wtTauMRCs{:,2});

[wtTauSort, tauInd] = sortrows(wtTauStats,1);
wtTauNameSort = wtTauNames(tauInd,:);
wtTauTraceSort = wtTauMRCTrace(tauInd,:);