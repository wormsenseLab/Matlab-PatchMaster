%metaDataConvert.m

function tracePicks = metaDataConvert(tracePicks)

if size(tracePicks, 2) == 3 %for trace list w/ cell, series, trace
    for i = 1:size(tracePicks,1)
        tracePicks{i,4} = str2num([tracePicks{i,3}]);
    end
    
    tracePicks = tracePicks(:, [1 2 4]);
    
elseif size(tracePicks,2) == 4
    for i = 1:size(tracePicks,1) %for calibTrace list w/ cell, series, trace, size
        tracePicks{i,5} = str2num([tracePicks{i,3}]);
        tracePicks{i,6} = str2num([tracePicks{i,4}]);
    end
    
    tracePicks = tracePicks(:, [1 2 5 6]);
end
    

end