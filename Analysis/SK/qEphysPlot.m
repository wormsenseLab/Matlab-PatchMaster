% qEphysPlot.m
%
% qEphysPlot (ephysData, cellName, protocol, channel, sweep)
% 
% Assumes channel 1 is current, plots pA instead of A

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

switch channel
    case 1
        plot(ephysData.(cellName).data{channel,protocol}(:,sweep)*1E12);
    otherwise
        plot(ephysData.(cellName).data{channel,protocol}(:,sweep));
end

end