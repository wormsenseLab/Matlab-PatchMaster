% FilterRecordings.m
% 
% 
% Created by Sammy Katta on 4 Jan 2018.

function filteredCells = FilterRecordings(ephysData,ephysMetadata,varargin)

p = inputParser;
% metadata spreadsheet (optional or required? only necessary if using the filtering parameters)
p.addRequired('ephysData');
p.addRequired('ephysMetadata'); 

p.addOptional('allCells', cell(0));

% If adding new filter parameters, beware the name must match the column
% name in RecordingDatabase.xlsx. Capitalization and spaces don't matter,
% but it must also be a valid struct fieldname.
% TODO: adjust the regexprep (~line 48) on metadata to get rid of special 
% characters (parentheses, symbols, greek letters) for full ability to
% match and filter by any column in the database.
p.addParameter('strain', cell(0));
p.addParameter('externalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('internalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('cellType', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('wormPrep', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));

p.parse(ephysData, ephysMetadata, varargin{:});

allCells = p.Results.allCells;

% doing it this way makes it easier to add new parameters straight into
% inputParser, rather than hardcoding parameter names.
filterParams = rmfield(p.Results,{'ephysData','ephysMetadata','allCells'});
filterHeaders = fieldnames(filterParams);

if isempty(allCells)
   allCells =  fieldnames(ephysData);
end

% find metadata rows for the cells specified by allCells to pull out the
% relevant parameters used for filtering.
allCellInd = cellfun(@(x) strcmp(ephysMetadata(:,1),x),allCells,'un',0)';
allCellInd = logical(sum([allCellInd{:}],2));


% find columns in metadata spreadsheet whose headers match the filter
% parameter names (case- and whitespace-insensitive).

useParams = [];

% only look at filter params that the user has specified
for i = 1:length(filterHeaders)  
    if ~isempty(filterParams.(filterHeaders{i}))
        useParams = [useParams i];
    end
end

paramInd = [];
subMetadata = cell(length(allCells),length(useParams));

for i=1:length(useParams)
    paramInd(i) = find(strcmpi(regexprep(ephysMetadata(1,:), '\s+', ''),filterHeaders{useParams(i)}));
    % convert numbers and NaNs (from empty slots) to strings to allow use 
    % of ismember for matching
    subMetadata(:,i) = cellfun(@num2str,ephysMetadata(allCellInd,paramInd(i)),'un',0);
end

% subMetadata = cellfun(@num2str,subMetadata,'un',0);
for i = 1:length(useParams)
    try isParam(:,i) = ismember(subMetadata(:,i),filterParams.(filterHeaders{useParams(i)}));
    catch
        fprintf('Could not match %s in metadata',char(filterParams.(filterHeaders{i})'));
        continue
    end
end

matchingRecs = sum(isParam,2)==length(useParams);

%NEXT: if param is a num, == it. 
%LATER: consider cases where you want to look at recordings within a range
%for a given param, e.g., 50<Rs<200MOhm. Might be better off doing that in
%Excel...or create a second filtering fxn.

filteredCells=allCells(matchingRecs);

end