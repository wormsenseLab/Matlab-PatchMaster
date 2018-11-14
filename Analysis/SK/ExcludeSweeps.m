% ExcludeSweeps.m
%
% selectedSweeps = ExcludeSweeps(ephysData, protList, varargin)
% 
% 
%TODO: inputparser param for electrode data vs. calib (plotting properly,
%plotting the right channel [either chan3 or search for 'mV'?)
%TODO: Build another quick GUI to look through tree of sweeps (maybe
%directory type? or 3 listboxes, first w fieldnames, second w series
%numbers that populate when you select a fieldname, third w channel
%number/channeltype/channelunit? single click? enter to plot?

function selectedSweeps = ExcludeSweeps(ephysData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('protList', @(x) iscell(x) && ~isempty(x) && ischar(x{1}));

p.addOptional('allCells', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}))

p.addParameter('channel', 1, @(x) isnumeric(x)); %1 for current, 2 for stim command, 3 for PD signal
p.addParameter('matchType', 'full', @(x) ischar(x));
p.addParameter('maxCols',4,@(x) isnumeric(x));
p.addParameter('addToList', cell(0), @(x) iscell(x));
p.addParameter('overwriteExisting',0,@(x) islogical(x) || isnumeric(x) && ismember(x,[0 1]));

p.parse(ephysData, protList, varargin{:});

allCells = p.Results.allCells;
channel = p.Results.channel;
matchType = p.Results.matchType;
maxCols = p.Results.maxCols;
existingSweeps = p.Results.addToList;
overwriteFlag = logical(p.Results.overwriteExisting);

if isempty(allCells)
    allCells = fieldnames(ephysData);
end

if isempty(existingSweeps)
    overwriteFlag = 1; %no reason to worry about overwriting prev. list if it's empty
end

maxPlots = 12;

protLoc = cell(length(allCells),1);
totSeries = 0;
selectedSweeps = cell(1,3);
if ~exist('matchType','var')
    matchType = 'full';
end

iCell = 1;
goBack = 0;
keepSweeps = cell(0);

while iCell <= length(allCells)
    protLoc{iCell} = matchProts(ephysData, allCells{iCell}, protList, 'matchType', matchType);
    thisCell = allCells{iCell};
    
    iSeries = 1;
    while iSeries <= length(protLoc{iCell})
        if goBack && lastSeries 
            iSeries = length(protLoc{iCell});
        end
        
        thisSeries = protLoc{iCell}(iSeries);
        data = ephysData.(allCells{iCell}).data{channel,protLoc{iCell}(iSeries)};
        dataType = ephysData.(allCells{iCell}).dataunit{channel,protLoc{iCell}(iSeries)};
        sf = ephysData.(allCells{iCell}).samplingFreq{protLoc{iCell}(iSeries)}/1000;
        protName = sprintf('%d: %s', iSeries,...
            ephysData.(thisCell).protocols{thisSeries});
        nSweeps = size(data,2);
        sweeps = 1:nSweeps;
        
        if ~overwriteFlag 
            %(if there were an existing list and overwriteFlag was true, it
            % would only overwrite those sweeps that overlapped but not everything)
           
%NEXT: check against oldSweeps list and increment iSeries? this has to work
%properly with goBack (i.e., repeat whatever the last GUI instance asked
%for, go one more previous or go one more next).
%Do this by reading in keepSweeps(oldSweepNums)=true if cellID/series
%matches and overwriteFlag=true. If not overwriteFlag, skip the series and
%proceed in the same direction as user had indicated (just don't reset goBack),
%and concatenate to selectedSweeps at the end, and sortrows there.

        end
        % Subtract the leak, and add it as dim 3 behind the raw trace for
        % easy passing to the GUI
        % use baseLength from first sweep of stimTree if available
        try baseLength = ephysData.(allCells{iCell}).stimTree{protLoc{iCell}(iSeries)}{3,3}.seDuration * 1e3 -1;
        catch
            baseLength = 30;
        end
        [leakSubtract, leakSize] = SubtractLeak(data,sf, 'BaseLength', baseLength);
        data(:,:,2) = leakSubtract;
        
        % Run the GUI, with a maximum number of plots per page (especially
        % useful for series with many reps). Can also adjust maximum
        % number of columns.

        iPage = 1;
        while iPage <= nSweeps
            
            if goBack == 1 && lastPage % if user hit previous on first page of new series, go to last page (of previous series)
                iPage = (floor(nSweeps/maxPlots))*maxPlots+1;
            end
            
            if nSweeps >= iPage+maxPlots-1
                pageSweeps = iPage:iPage+maxPlots-1;
            else
                pageSweeps = iPage:nSweeps;
            end
            
            if nSweeps < maxCols
                nCols = nSweeps;
            else
                nCols = maxCols;
            end
            
            pageNo = [pageSweeps(1) pageSweeps(end) nSweeps];
            try keepPageSweeps = keepSweeps{iCell}{iSeries}(pageSweeps);  
            catch % if keepSweeps hasn't been initialized here, initialize as true
                keepSweeps{iCell}{iSeries} = true(size(sweeps)); % assume all sweeps are kept unless user excludes them
                keepPageSweeps = keepSweeps{iCell}{iSeries}(pageSweeps);
            end
            
            % Run the GUI for the current series
            %TODO: Get Esc key to pass out an error to catch (currently, error
            %seems to be passed to uiwait instead of out to ExcludeSweeps).
            %TODO: Get -1,0,+1 output for next/previous button and use to
            %modify iSeries. Okay to clear previous selection, or do we need
            %to replay excluded traces?
            [keepPageSweeps, goBack] = selectSweepsGUI(...
                data(:,pageSweeps,:),dataType,channel,leakSize(pageSweeps),sf,thisCell,protName,nCols,pageNo,keepPageSweeps);
            
            % indexing sweeps allows GUI to remember which sweeps were
            % selected between pages in a recording, but will reset across
            % series.
            keepSweeps{iCell}{iSeries}(pageSweeps) = keepPageSweeps; % 

            if goBack && iPage == 1 % if user hit previous on first page of new series, go to previous series
                lastPage = 1;
                break
            elseif goBack % if user hit previous any other time, go back a page
                iPage = iPage - maxPlots;
                lastPage = 0;
            else % if user hit enter/next, continue onward and save the selected sweeps
                iPage = iPage + maxPlots;
                lastPage = 0;
            end
        end
        
        % Format the list of sweeps into a text string for the Excel sheet
        keepSweepNums{iCell}{iSeries} = sweeps(keepSweeps{iCell}{iSeries});

        if iSeries == 1 && goBack
            iCell = max([iCell - 1, 1]); %either subtract one or stay if this is first series
            lastSeries = 1;
            break
        elseif goBack
            iSeries = iSeries - 1;
            lastSeries = 0;
        else
            iSeries = iSeries+1;
            lastSeries = 0;
        end
    end
    
    
    if iCell == 1 && goBack
        fprintf('First series, please select traces to exclude.');
        continue
    elseif goBack
        iCell = iCell - 1;
    else
        iCell = iCell + 1;
    end
    
    
end

keeps = ~cellfun(@isempty,keepSweepNums);

keepCells = allCells(keeps);
keepProts = protLoc(keeps);
keepSweepNums = keepSweepNums(keeps);


% Write out the kept sweep numbers into a single cell array for exporting
% to Excel.
totSeries = 0;
for iCell = 1:length(keepCells)
    for iSeries = 1:length(keepProts{iCell})
        totSeries = totSeries + 1;
        sweeps = keepSweepNums{iCell}{iSeries};
        
        if ~isempty(sweeps)
            
            sweepsTxt = num2str(sweeps,'%g,'); % add commas within string
            sweepsTxt = sweepsTxt(1:end-1); % trim last comma
            sweepsTxt = horzcat('''',sweepsTxt); % prepend ' to force xls text format
            
            selectedSweeps{totSeries,1}=keepCells{iCell};
            selectedSweeps{totSeries,2}=keepProts{iCell}(iSeries);
            selectedSweeps{totSeries,3}=sweepsTxt;

        end
    end
end

keepRows = cellfun(@(x) ~isempty(x), selectedSweeps(:,3));
selectedSweeps = selectedSweeps(keepRows,:);

% Ask where to save and write out the .xls file)
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save sweep list to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName, selectedSweeps);


end