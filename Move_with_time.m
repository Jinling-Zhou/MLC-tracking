% Move MLC with time.
% Jinling Zhou, 4/3/2023.

clc;
clear all;
close all;

% Read leaf positions and boundaries
path = 'RP.QA303005.test.dcm';
DicomInfo = dicominfo(path);
LeafJawPositions = DicomInfo.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions;
LeafPositionBoundaries = DicomInfo.BeamSequence.Item_1.BeamLimitingDeviceSequence.Item_3.LeafPositionBoundaries;

Leaf_A_o = LeafJawPositions(1:60);
Leaf_B_o = LeafJawPositions(61:120);
for i = 1: length(LeafPositionBoundaries)-1
center_value(i) = (LeafPositionBoundaries(i)+LeafPositionBoundaries(i+1))/2;
end

% Read the movement data
Movement = readmatrix("sine wave_A2T3.csv");

t= Movement(:,1);
Shift_t= Movement(:,2);

% Move MLC as described by Movement 
for i = 1:length(t)
    Leaf_D_Shift = Shift_t(i)* cos(-pi/2)*10;
    Edge_D_Shift = Shift_t(i)* sin(-pi/2)*10;
    [A_f, B_f] = singlemove(Leaf_A_o, Leaf_B_o,LeafPositionBoundaries,Leaf_D_Shift,Edge_D_Shift);
    A_final(i,:) = A_f;
    B_final(i,:) = B_f;
    center_value_t(i,:) = center_value;
    for j = 1: length(center_value)
        t_space(i,j) = t(i);
    end    
end
 plot3(A_final,center_value_t,t_space,'*m',B_final,center_value_t,t_space,'*m');
%  axis equal;
daspect([10 10 1]);
 hold off

% Single movement of bank A and bank B
function [Leaf_A_final, Leaf_B_final] = singlemove(Leaf_A_original, Leaf_B_original, LeafPositionBoundaries, LeafDirectionShift, EdgeDirectionShift)
% Move along leaf direction:
Leaf_A_final = Leaf_A_original + LeafDirectionShift;
Leaf_B_final = Leaf_B_original + LeafDirectionShift;
Leaf_A_move_Parallel = Leaf_A_original + LeafDirectionShift;
Leaf_B_move_Parallel = Leaf_B_original + LeafDirectionShift;
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
       if i+j <=60 && i+j >=1
        Leaf_A_final(i+j) = Leaf_A_move_Parallel(i);
        Leaf_B_final(i+j) = Leaf_B_move_Parallel(i);
       end 
end       
end
