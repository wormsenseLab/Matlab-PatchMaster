% peakRatios.m
% 
% This function calculates the ratio between two peaks either in the same
% trace or different traces, given the location timepoint.
% 
% USAGE:
%   test = peakRatios(p1Time, p2Time, trace)
%   test = peakRatios(p1Time, p2Time, trace1, trace2)
% 


function peakRatios(seriesStimuli, varargin)

% If second peak is in a separate trace, read that in instead
% TODO: rewrite this w inputparser, if no trace2 then copy trace1 to trace2

p = inputParser;
p.addRequired('seriesStimuli');
p.addRequired('cellName');
p.addRequired('protName');

p.addParameter('MatchType', 'full', @(x) validateattributes(x,{'char'},{'nonempty'}));
p.parse(ephysData, cellName, protName, varargin{:});

matchType = p.Results.MatchType;
    
end
