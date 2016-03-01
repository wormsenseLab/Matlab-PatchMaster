% ExcludeSweeps.m
%
%

% add matchtype to pass through to matchProts?
% pass channel too?
function selectedSweeps = ExcludeSweeps(ephysData, allCells, protList)

protLoc = cell(length(allCells),1);
totSeries = 0;

for iCell = 1:length(allCells)
    protLoc{iCell} = matchProts(ephysData, allCells{iCell}, protList);
    cellName = allCells{iCell};
    
    for iSeries = 1:length(protLoc{iCell})
        data = ephysData.(allCells{iCell}).data{1,protLoc{iCell}(iSeries)};
        sf = ephysData.(allCells{iCell}).samplingFreq{protLoc{iCell}(iSeries)}/1000;
        nSweeps = size(data,2);
        sweeps = 1:nSweeps;
        
        [leakSubtract, leakSize] = SubtractLeak(data,sf);
        data(:,:,2) = leakSubtract;
        try keepSweeps = selectSweepsGUI(data,leakSize);
        catch
            fprintf('Exited on %s series %d',cellName,protLoc{iCell});
            return;
        end
        
        
        % Format the list of sweeps into a text string for the Excel sheet
        sweeps = sweeps(keepSweeps);
        
        if ~isempty(sweeps)
            totSeries = totSeries+1;
            
            sweepsTxt = num2str(sweeps,'%g,'); % add commas within string
            sweepsTxt = sweepsTxt(1:end-1); % trim last comma
            sweepsTxt = horzcat('''',sweepsTxt); % prepend ' to force xls text format
            % (necessary when only 1 sweep left)
            
            % Set up the cell name, series number, and sweep numbers in the
            % right format for outputting to the Excel sheet
            selectedSweeps{totSeries,1}=cellName;
            selectedSweeps{totSeries,2}=protLoc{iCell}(iSeries);
            selectedSweeps{totSeries,3}=sweepsTxt;
        end
        
    end
    
    
end

% Ask where to save and write out the .xls file)
[filename, pathname] = uiputfile(...
    {'*.xls;*.xlsx', 'Excel files';
    '*.*', 'All files'}, ...
    'Save sweep list to .xls file:');
fName = fullfile(pathname,filename);
xlswrite(fName, selectedSweeps);


end