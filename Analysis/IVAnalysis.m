 % IVAnalysis.m
%
% IVAnalysis takes the data from on-cell and whole-cell IV steps, subtracts
% on-cell from whole-cell as capacitance correction, and takes the mean of
% three subsequent subtracted

function capCorrIV = IVAnalysis(ephysData, varargin)

nGroups = nargin-1;
nReps = 3;

capCorrIV = cell(1,nGroups);

if nGroups < 1
    error('NoGroups','Please provide names of recordings in a columnar cell array.')
end


for iGroup = 1:nGroups
    groupCells = varargin{iGroup};
    
    for iCell = 1:length(groupCells)
        cellName = groupCells{iCell}; %TODO: maybe split into project name and cell numbers when feeding input
        
        % Find locations of protocols ending in 'IVq'.
        protName = 'IVq';
        protLoc = matchProts(ephysData,cellName,protName,'MatchType','last');
        
        % Only run analysis if this recording has IVq protocols
        if protLoc
            % Define the positions of the IVq series for analysis, out of
            % the set of recorded IVqs for that cell. (i.e, if protWC = 4,
            % it's the 4th IVq series for the cell, regardless of OC vs. WC
            % or whether other series exist in between).
            % Currently set by parent script.
            protOC = ephysData.(cellName).protOC;
            protWC = ephysData.(cellName).protWC;
            
            for i = 1:length(protLoc)
                ivq{i} = ephysData.(cellName).data{1,protLoc(i)};
            end
            % Take all the iv steps gathered for the given cell
            ivSteps{iCell} = ivq;
            
            if protOC ~= 0 % only do cap correct if OC_ivq exists
                % Pull out the desired two OC and WC iv steps and save by cell name
                % to make it more user friendly when later separating/replotting
                % particular cells of interest. Turn cells into matrix, then split
                % matrix up so you can later take the mean of the three sets of
                % capacitance-corrected iv steps.
                currSteps = ivSteps{iCell}(protOC:protOC+nReps-1);
                protSize = size(ivSteps{iCell}{protOC});
                eval(sprintf('%s_OC = reshape(cell2mat(currSteps),[protSize nReps]);',cellName))
                currSteps = ivSteps{iCell}(protWC:protWC+nReps-1);
                eval(sprintf('%s_WC = reshape(cell2mat(currSteps),[protSize nReps]);',cellName))
                
                
                % Capacitance correction by subtracting the on-cell from the
                % whole-cell (both as the mean of three technical replicates)
                eval(sprintf('capCorrCells.(cellName).capCorrIV = mean(%s_WC,3)-mean(%s_OC,3);',cellName,cellName))
            else
                protSize = size(ivSteps{iCell}{protWC});
                currSteps = ivSteps{iCell}(protWC:protWC+nReps-1);
                eval(sprintf('%s_WC = reshape(cell2mat(currSteps),[protSize nReps]);',cellName))                
                eval(sprintf('%s_OC = zeros(2,2,2);',cellName))
                
                eval(sprintf('capCorrCells.(cellName).capUnCorrIV = mean(%s_WC,3);',cellName))
            end
            
        end
        
        eval(sprintf('capCorrCells.(cellName).OCmean = mean(%s_OC,3);', cellName))
        eval(sprintf('capCorrCells.(cellName).WCmean = mean(%s_WC,3);', cellName))
        
    end
    
    capCorrIV{iGroup} = capCorrCells;
    clear capCorrCells;
end



end
