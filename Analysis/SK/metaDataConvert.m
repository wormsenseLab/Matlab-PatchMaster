%metaDataConvert.m

function tracePicks = metaDataConvert(tracePicks)

for i = 1:size(tracePicks,1)
    tracePicks{i,4} = str2num([tracePicks{i,3}]);
end

tracePicks = tracePicks(:, [1 2 4]);

end