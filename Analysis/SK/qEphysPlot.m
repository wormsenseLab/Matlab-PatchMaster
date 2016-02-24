% qEphysPlot.m
% 
% qEphysPlot (ephysData, cellName, protocol, channel, sweep)

function qEphysPlot (ephysData, cellName, protocol, varargin)

p = inputParser;
p.addRequired('ephysData');
p.addRequired('cellName');
p.addRequired('protocol');

p.addOptional('channel', 1);
p.addOptional('sweep',1);

p.parse(ephysData,cellName, protocol, varargin{:});
channel = p.Results.channel;
sweep = p.Results.sweep;

plot(ephysData.(cellName).data{channel,protocol}(:,sweep));

end