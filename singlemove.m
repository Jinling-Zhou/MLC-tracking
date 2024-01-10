% This program moves MLC in any direction, for one single movement only.
% Jinling Zhou, 4/3/2023.
clc;
clear all;
close all;

% Read the leaf positions and boundary locations
path = 'RP.QA303005.test.dcm';
DicomInfo = dicominfo(path);
LeafJawPositions = DicomInfo.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions;
LeafPositionBoundaries = DicomInfo.BeamSequence.Item_1.BeamLimitingDeviceSequence.Item_3.LeafPositionBoundaries;

% Separate the MLC into bank A and bank B
Leaf_A_original = LeafJawPositions(1:60);
Leaf_B_original = LeafJawPositions(61:120);

% Movement vector
LeafDirectionShift = 10 * cos(pi/2);
EdgeDirectionShift = 10 * sin(pi/2);

% Visualization
% f1 = figure;
% plot(LeafJawPositions);
% f2 = figure;
% plot(LeafPositionBoundaries);

% trace.txt was created to debug
% fileID = fopen('trace.txt', 'a');

% Move along leaf direction:
Leaf_A_move_Parallel = Leaf_A_original + LeafDirectionShift;
Leaf_B_move_Parallel = Leaf_B_original + LeafDirectionShift;
Leaf_A_final = Leaf_A_move_Parallel;
Leaf_B_final = Leaf_B_move_Parallel;
% Move along the edge direction:
for i = 1: length(LeafPositionBoundaries)-1
% Move perpendicular to leaf direciton, j is the number of edge locations
% found within EdgeDirectionShift
        center_value(i) = (LeafPositionBoundaries(i)+LeafPositionBoundaries(i+1))/2;
        a = center_value (i);
        b = a + EdgeDirectionShift;   
        if EdgeDirectionShift >= 0
        j = sum(LeafPositionBoundaries > a & LeafPositionBoundaries < b);
        else
        j = -sum(LeafPositionBoundaries > b & LeafPositionBoundaries < a);
        end
%         fprintf(fileID,'i=%d, j=%d, i+j=%d\n', i, j, i+j);
       if i+j <=60 && i+j >=1
        Leaf_A_final(i+j) = Leaf_A_move_Parallel(i);
        Leaf_B_final(i+j) = Leaf_B_move_Parallel(i);
       end 
 end       
% fclose(fileID);

% plot
figure;
x = center_value(1:60);
plot(x,Leaf_A_original,'og',x,Leaf_B_original,'og',x,Leaf_A_final,'*m',x,Leaf_B_final,'*m');
% axis([-200,200, -200,200]);
axis equal;
