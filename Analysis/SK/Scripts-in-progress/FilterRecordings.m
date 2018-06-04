% FilterRecordings.m
% 
% This function filters the list of recording names to be analyzed based on
% metadata in the RecordingDatabase.xlsx format. 
% 
% USAGE: 
%   filteredCells = FilterRecordings(ephysData,ephysMetadata)
%   filteredCells = FilterRecordings(ephysData,ephysMetadata,allCells)
%   filteredCells = FilterRecordings(ephysData,ephysMetadata,'internal',{'IC2'})
% 
% INPUTS:
%   ephysData           struct      Struct containing data created by the
%                                   function ImportPatchData.
% 
%   ephysMetadata       cell        Cell array version of metadata records,
%                                   created by using ImportMetaData() on
%                                   RecordingDatabase.xlsx.
% 
% OPTIONAL INPUTS:
%   allCells            cell        Cell array of strings, containing a
%                                   list of recording names. Useful if
%                                   further narrowing down an
%                                   already-filtered list of recordings.
% 
% FILTER PARAMETERS:
%   Parameter-value pairs can be used to specify which columns of the
%   recording database to filter by, such as strain name or solutions used.
%   The parameter name references the column (e.g., 'strain'), and the 
%   value should be a cell array of strings (e.g., {'TU2769'} or
%   {'TU2769','TU253'} if you are looking to match values. 
% 
%   To filter using numerical ranges, the value can be specified either as
%   a lower/upper bound in a 1x2 vector, e.g., [40 80]. Alternatively, you
%   may use a string specifying less than/equal to/greater than with a
%   single number, e.g., '<250' or '=45'. By default, the function includes
%   the specified number (i.e., '<250' is actually <=250). 
%   When using equals, check that your spreadsheet is not rounding numbers
%   in the display. You may also provide a single number (a 1x1 double) for
%   equals, which will be matched as a string.
%   
%   Current string parameter options: 'strain', 'internalSolution',
%   'externalSolution', 'cellType', 'wormPrep'.
% 
%   Current numerical parameter options: 'stimFilterFrequencykHz', 
%   'cellStimDistUm', 'RsM', 'TCultC', 'roomTempC'.
% 
% 
% 
%   The function can be expanded to include other string-based filters by
%   simply adding an addParameter line to the inputParser section at the
%   beginning. The parameter name must match the name of the column in the
%   metadata sheet (case-insensitive, ignoring special characters if the
%   ignoreSpecial parameter is set to true.
% 
% OTHER PARAMETERS:
%   ignoreSpecial       true/false  If set to true, the function will
%                                   ignore all non-alphanumeric characters
%                                   (no symbols, Greek letters) when 
%                                   matching filter parameters to columns.
% 
% 
% 
% OUTPUTS:
%   filteredCells       cell array  Cell array of strings with names of
%                                   recordings that match all the given
%                                   filter parameters.
% 
%   filteredMetadata    cell array  The subset of the metadata spreadsheet
%                                   corresponding to filteredCells.
% 
% Created by Sammy Katta on 4 Jan 2018.

function [filteredCells, filteredMetadata] = FilterRecordings(ephysData,ephysMetadata,varargin)

% PARSE INPUTS
p = inputParser;
% metadata spreadsheet (optional or required? only necessary if using the filtering parameters)
p.addRequired('ephysData');
p.addRequired('ephysMetadata'); 

p.addOptional('allCells', cell(0));
p.addParameter('ignoreSpecial', true);


% If adding new filter parameters, beware the name must match the column
% name in RecordingDatabase.xlsx. Capitalization and spaces don't matter,
% but it must also be a valid struct fieldname.
% TODO: adjust the regexprep (~line 48) on metadata to get rid of special 
% characters (parentheses, symbols, greek letters) for full ability to
% match and filter by any column in the database.
p.addParameter('strain', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('externalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('internalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('cellType', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('wormPrep', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('stimLocation', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1})); % anterior or posterior
p.addParameter('cellStimDistUm',cell(0)); % distance of stim from cell body, um
p.addParameter('stimFilterFrequencykHz',cell(0)); % external stimulus filter frequency, kHz
p.addParameter('roomTempC',cell(0)); % room temperature, degC
p.addParameter('TCultC',cell(0)); % cultivation temperature, degC
p.addParameter('RsM',cell(0)); % series resistance, in MOhms


p.parse(ephysData, ephysMetadata, varargin{:});

allCells = p.Results.allCells;
stripFlag = logical(p.Results.ignoreSpecial);

% doing it this way makes it easier to add new parameters straight into
% inputParser, rather than hardcoding parameter names.
filterParams = rmfield(p.Results,{'ephysData','ephysMetadata','allCells','ignoreSpecial'});
filterHeaders = fieldnames(filterParams);

if isempty(allCells)
   allCells =  fieldnames(ephysData);
end


% FILTERING CODE STARTS HERE
% find metadata rows for the cells specified by allCells to pull out the
% relevant parameters used for filtering.
allCellInd = cellfun(@(x) find(strcmp(ephysMetadata(:,1),x)),allCells,'un',0)';
allCellInd = [allCellInd{:}];


% find columns in metadata spreadsheet whose headers match the filter
% parameter names (case- and whitespace-insensitive).

useParams = [];
paramType = cell(0);

% only look at filter params that the user has specified
for i = 1:length(filterHeaders) 
    thisParam = filterParams.(filterHeaders{i});
    if ~isempty(thisParam)        
        % determine whether user has given a string, or specified a
        % numerical range in either string or vector format.
        switch class(thisParam)
            case 'cell'
                isRange = regexp(thisParam,'[<=>]');
                isRange = [isRange{:}];
                
                if isRange
                    paramType = horzcat(paramType,'textrange');
                else
                    paramType = horzcat(paramType,'string');
                end
                
            case 'char'
                filterParams.(filterHeaders{i}) = {filterParams.(filterHeaders{i})};
                thisParam = {thisParam};

                isRange = regexp(thisParam,'[<=>]');
                isRange = [isRange{:}];
                
                if isRange
                    paramType = horzcat(paramType,'textrange');
                else
                    paramType = horzcat(paramType,'string');
                end
                
            case 'double'
                if size(thisParam,2) == 2
                    paramType = horzcat(paramType,'range');
                elseif size(thisParam,1) == 1
                    filterParams.(filterHeaders{i}) = {num2str(filterParams.(filterHeaders{i}))};
                    paramType = horzcat(paramType,'string');
                else
                    fprintf('Skipped %s filter - please format numerical filters as specified in function help', filterHeaders(i));
                    continue
                end

        end
        
        useParams = [useParams i]; % mark this parameter as user-specified
    end
end

paramInd = [];
subMetadata = cell(length(allCells),length(useParams));


% Pull out metadata columns specified for filtering, using case-insensitive
% string comparison. Ignore only whitespace if stripFlag (ignoreSpecial) is
% is not set.

%TEST: parameter columns with special characters
for i=1:length(useParams)
    if stripFlag
        paramInd(i) = find(strncmpi(regexprep(ephysMetadata(1,:), '[^a-zA-Z0-9]', ''),...
            regexprep(filterHeaders{useParams(i)}, '[^a-zA-Z0-9]', ''),length(filterHeaders{useParams(i)})));
        
    else
        paramInd(i) = find(strncmpi(regexprep(ephysMetadata(1,:), '\s+', ''),filterHeaders{useParams(i)},length(filterHeaders{useParams(i)})));
    end
    % if the given parameter is a string to be matched, make the sure the
    % metadata is a string. If it's a range (whether specified numerically
    % or by text), make sure the metadata is double or NaN.
    switch paramType{i}
        case 'string'
            subMetadata(:,i) = cellfun(@num2str,ephysMetadata(allCellInd,paramInd(i)),'un',0);
            try isParam(:,i) = ismember(lower(subMetadata(:,i)),lower(filterParams.(filterHeaders{useParams(i)})));
            catch
                fprintf('Could not match %s in metadata',char(filterParams.(filterHeaders{i})'));
                continue
            end
            
            
        case 'range'
            subMetadata(:,i) = cellfun(@double,ephysMetadata(allCellInd,paramInd(i)),'un',0);
            try isParam(:,i) = filterParams.(filterHeaders{useParams(i)})(1) <= [subMetadata{:,i}] ...
                    & filterParams.(filterHeaders{useParams(i)})(2) >= [subMetadata{:,i}];
            catch
                fprintf('Error matching range for %s in metadata\n',char(filterParams.(filterHeaders{i})'));
            end
            
        % if mathematical range was given as text, do the appropriate
        % comparison. 
        case 'textrange'
            subMetadata(:,i) = cellfun(@double,ephysMetadata(allCellInd,paramInd(i)),'un',0);
            mathString = char(filterParams.(filterHeaders{useParams(i)}));
            
            lessMatch = ~isempty(regexp(mathString, '<'));
            equalMatch = ~isempty(regexp(mathString, '='));
            greaterMatch = ~isempty(regexp(mathString,'>'));
            
            mathParam = regexp(mathString,'[<=>]','split');
            mathParam = str2num(mathParam{2});
            
            if equalMatch && ~lessMatch && ~greaterMatch
                try isParam(:,i) = [subMetadata{:,i}] == mathParam;
                catch
                    fprintf('Error matching range for %s in metadata\n',char(filterParams.(filterHeaders{i})'));
                end
            elseif lessMatch && ~greaterMatch
                try isParam(:,i) = [subMetadata{:,i}] <= mathParam;
                catch
                    fprintf('Error matching range for %s in metadata\n',char(filterParams.(filterHeaders{i})'));
                end
            elseif greaterMatch
                try isParam(:,i) = [subMetadata{:,i}] >= mathParam;
                catch
                    fprintf('Error matching range for %s in metadata\n',char(filterParams.(filterHeaders{i})'));
                end
            else
                fprintf('Error matching range for %s in metadata.\nPlease refer to function help for proper formatting of numerical range specification.',char(filterParams.(filterHeaders{i})'))
                
            end
    end
    
    
    
end

% Find which recordings have parameters matching all of the individual
% specified parameters (case-insensitive for strings).

matchingRecs = sum(isParam,2)==length(useParams);

filteredCells=allCells(matchingRecs);

% Optional output: the complete metadata subset for just those matching cells.
filteredMetadata = ephysMetadata(allCellInd,:);
filteredMetadata = filteredMetadata(matchingRecs,:);

end