% matchProts.m
%
% Find the locations/indices of series matching the desired protocol name.
% matchProts is case-insensitive.
%
% USAGE:
%   protLoc = matchProts(ephysData, cellName, protName)
%   protLoc = matchProts(ephysData, cellName, protName, 'MatchType', 'first')
%
% INPUTS:
%   ephysData       struct      Imported data from ImportPatchData.
%
%   cellName        char        Name of the recording to look at (must be a
%                               field within ephysData).
%
%   protName        char/cell   String or cell array of strings with the
%                               names of protocols to be matched.
%
% OUTPUTS:
%   protLoc         double      Vector of the indices of all series that
%                               match the given set of protocol names, in
%                               numerical order.
%
% OPTIONAL PARAMETERS:
%   'MatchType;     char        String specifying how to match the protocol
%                               name. Defaults to 'full'.
%                                   'full' finds series where the name
%                                       exactly matches protName
%                                   'first' finds series where the start of
%                                       the series name matches protName
%                                   'last' finds series where the end of
%                                       the series name matches protName
%
% Created by Sammy Katta on 23 February 2016.

function protLoc = matchProts(ephysData, cellName, protName, varargin)

% Parse inputs and allow proper input of parameter/value pairs.
p = inputParser;
p.addRequired('ephysData');
p.addRequired('cellName');
p.addRequired('protNames');

p.addParameter('MatchType', 'full', @(x) validateattributes(x,{'char'},{'nonempty'}));
p.parse(ephysData, cellName, protName, varargin{:});

matchType = p.Results.MatchType;

% If only one protocol name was given, make it a cell.
if ischar(protName) && ~iscell(protName)
    protName = {protName};
end

% Find matches for the protocol name in the list of series names for the
% given recording, based on the user's specification.
protLoc = cell(length(protName),1);
for i = 1:length(protName)
    switch matchType
        case 'full' % match the string exactly
            prots = ephysData.(cellName).protocols;
            protLoc{i} = find(strcmpi(protName{i},prots));
        case 'first' % match protName to the beginning of the series name
            prots = ephysData.(cellName).protocols;
            protLoc{i} = find(strncmpi(protName{i},prots,length(protName)));
        case 'last' % match protName to the end of the series name
            flippedProts = cellfun(@fliplr, ephysData.(cellName).protocols, ...
                'UniformOutput', false);
            protLoc{i} = find(strncmpi(fliplr(protName{i}),flippedProts,length(protName)));
    end
    
end

% Combine the list of all found series indices into one double vector,
% sorted numerically. Indices that matched with multiple protocol names are
% only included once.
protLoc = unique(sort([protLoc{:}]));

end