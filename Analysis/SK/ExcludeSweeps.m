% ExcludeSweeps.m
%
% selectedSweeps = ExcludeSweeps(ephysData, allCells, channel, protList, matchType)
% 
% 
%TODO: inputparser for matchtype
%TODO: inputparser param for electrode data vs. calib (plotting properly,
%plotting the right channel [either chan3 or search for 'mV'?)
%TODO: Build another quick GUI to look through tree of sweeps (maybe
%directory type? or 3 listboxes, first w fieldnames, second w series
%numbers that populate when you select a fieldname, third w channel
%number/channeltype/channelunit? single click? enter to plot?

function selectedSweeps = ExcludeSweeps(ephysData, allCells, channel, protList, matchType)

protLoc = cell(length(allCells),1);
totSeries = 0;
selectedSweeps = cell(1,3);
if ~exist('matchType','var')
    matchType = 'full';
end

for iCell = 1:length(allCells)
    protLoc{iCell} = matchProts(ephysData, allCells{iCell}, protList, 'matchType', matchType);
    thisCell = allCells{iCell};
    
    wSeries = 1;
    while wSeries <= length(protLoc{iCell})
        thisSeries = protLoc{iCell}(wSeries);
        data = ephysData.(allCells{iCell}).data{channel,protLoc{iCell}(wSeries)};
        dataType = ephysData.(allCells{iCell}).dataunit{channel,protLoc{iCell}(wSeries)};
        sf = ephysData.(allCells{iCell}).samplingFreq{protLoc{iCell}(wSeries)}/1000;
        protName = sprintf('%d: %s', wSeries,...
            ephysData.(thisCell).protocols{thisSeries});
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
        [keepSweeps, goBack] = selectSweepsGUI(data,dataType,channel,leakSize,sf,thisCell,protName);
        %         catch
        %             fprintf('Exited on %s series %d',cellName,protLoc{iCell});
        %             return;
        %         end
        
        
        % Format the list of sweeps into a text string for the Excel sheet
        sweeps = sweeps(keepSweeps);
        
        if ~isempty(sweeps)
            
            
%TODO: find whether an entry exists for the current series. If so,
%overwrite it. If not, make a new one (at the end?)
            %
%             if ~ismember([selectedSweeps{ismember(selectedSweeps(:,1),thisCell),2}],thisSeries)
                totSeries = totSeries+1;
%             else
%                 totSeries = find(ismember(selectedSweeps(:,1),thisCell) && ismember([selectedSweeps{:,2}],thisSeries));
%             end
%             
            sweepsTxt = num2str(sweeps,'%g,'); % add commas within string
            sweepsTxt = sweepsTxt(1:end-1); % trim last comma
            sweepsTxt = horzcat('''',sweepsTxt); % prepend ' to force xls text format
            % (necessary when only 1 sweep left)
            
            % Set up the cell name, series number, and sweep numbers in the
            % right format for outputting to the Excel sheet
            selectedSweeps{totSeries,1}=thisCell;
            selectedSweeps{totSeries,2}=thisSeries;
            selectedSweeps{totSeries,3}=sweepsTxt;
        end
        
        % If user has clicked previous button, display the previous series.
%TODO: Fix this to jump back within wSeries and totSeries and iCell.
%Luckily, you sillily used totSeries instead of just using end, so now you
%can make use of that. Hurray!
        if goBack && totSeries > 1
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