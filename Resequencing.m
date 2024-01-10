clc;
clear;
close all;

% Original data
load('Result_Perpendicular.mat');
givenMatrix = Result_Perpendicular;

% Get unique values of i, j, k
unique_i = unique(givenMatrix(:, 1));
unique_j = unique(givenMatrix(:, 2));
unique_k = unique(givenMatrix(:, 3));

% Preallocate the new matrix
newMatrix = zeros(size(givenMatrix));

m = 1;
% Loop through unique values of i
for i_idx = 1:numel(unique_i)
    current_i = unique_i(i_idx);
    
    % Loop through unique values of j
    for j_idx = 1:numel(unique_j)
        current_j = unique_j(j_idx);
        
        % Loop through unique values of k in reverse order
        for k_idx = numel(unique_k):-1:1
            current_k = unique_k(k_idx);
            
            % Find rows in the given matrix that match the current i, j, k
            idx = givenMatrix(:, 1) == current_i & givenMatrix(:, 2) == current_j & givenMatrix(:, 3) == current_k;
            
            % Assign the corresponding row from the given matrix to the new matrix
            newMatrix(m, :) = givenMatrix(idx, :);
            m = m+1;
        end
    end
end

