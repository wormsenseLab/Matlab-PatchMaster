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
    
    %TODO: change this to while loop to get Previous button to work (allows
    %index jumping)
    wSeries = 1;
    while wSeries <= length(protLoc{iCell})
        data = ephysData.(allCells{iCell}).data{1,protLoc{iCell}(wSeries)};
        sf = ephysData.(allCells{iCell}).samplingFreq{protLoc{iCell}(wSeries)}/1000;
        protName = sprintf('%d: %s', wSeries,...
            ephysData.(allCells{iCell}).protocols{protLoc{iCell}(wSeries)});
        nSweeps = size(data,2);
        sweeps = 1:nSweeps;
        
        % Subtract the leak, and add it as dim 3 behind the raw trace for
        % easy passing to the GUI
        [leakSubtract, leakSize] = SubtractLeak(data,sf);
        data(:,:,2) = leakSubtract;
        
        % Run the GUI for the current series
        %TODO: Get Esc key to pass out an error to catch (currently, error
        %seems to be passed to uiwait instead of out to ExcludeSweeps).
        %TODO: Get -1,0,+1 output for next/previous button and use to
        %modify iSeries. Okay to clear previous selection, or do we need
        %to replay excluded traces?
        [keepSweeps, goBack] = selectSweepsGUI(data,leakSize,cellName,protName);
%         catch
%             fprintf('Exited on %s series %d',cellName,protLoc{iCell});
%             return;
%         end
        
        
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
            selectedSweeps{totSeries,2}=protLoc{iCell}(wSeries);
            selectedSweeps{totSeries,3}=sweepsTxt;
        end
        
        if goBack && wSeries > 1
            wSeries = wSeries-1;
        elseif goBack && wSeries <= 1
            fprintf('On first series, can''t go back.');
        else
            wSeries = wSeries+1;
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