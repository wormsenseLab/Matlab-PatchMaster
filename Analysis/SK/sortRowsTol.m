% sortRowsTol.m
% 
% This functions sorts a matrix and finds unique rows given a tolerance
% value for each column.
% 
% USAGE:
% 
% 
% Created by Sammy Katta on 12 June 2018.
% Based on:
% https://stackoverflow.com/questions/19892549/return-unique-rows-with-a-tolerance-for-a-matrix

function [matA_sorted, sortIdx, matA_unique, startIdx, endIdx] = sortRowsTol(matA, tol, sortOrder)

[matA_sorted, sortIdx] = sortrows(matA);

for iCol = 1:length(sortOrder)
    
    %pick out values which have a diff greater than tolerance from their
    %neighbor after sorting. 
    isUnq(:,iCol) = [1; abs(diff(abs(matA_sorted(:,sortOrder(iCol)))))>tol(iCol)];
    isUnqEnd(:,iCol) = [abs(diff(abs(matA_sorted(:,sortOrder(iCol)))))>tol(iCol);1];
    
end

% These need to be indices in the sorted matrix to be useful.
startIdx = find(logical(any(isUnq,2)));
endIdx = find(logical(any(isUnqEnd,2)));

matA_unique = matA_sorted(startIdx,:); % result

end
