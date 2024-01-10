% This file moves MLC (a single move, in positive directions) and write the plan back to a DICOM file.
% for moving directions towards the leaf edges, can only shift to the right/positive/higher MLC number direction. 
% For fully functional single movement at any directions, check SingleMove.m.
% Jinling Zhou, 3/29/2023.
clc;
clear all;
close all;

% read the DICOM file.
path = 'RP.QA303005.test.dcm';
DicomInfo = dicominfo(path);
LeafJawPositions = DicomInfo.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions;
LeafPositionBoundaries = DicomInfo.BeamSequence.Item_1.BeamLimitingDeviceSequence.Item_3.LeafPositionBoundaries;
LeafDirectionShift = 50 * cos(pi/4);
EdgeDirectionShift = 50 * sin(pi/4);

% Visualization
% f1 = figure;
% plot(LeafJawPositions);
% f2 = figure;
% plot(LeafPositionBoundaries);

% Move along leaf direction:
Leaf_Positions_Parallel = LeafJawPositions + LeafDirectionShift;
for i = 1: length(LeafPositionBoundaries)-1
% Move perpendicular to leaf direciton, j is the number of edge locations
% found within EdgeDirectionShift
        a = LeafPositionBoundaries(i);
        b = LeafPositionBoundaries(i) + EdgeDirectionShift;
        center_value(i)=(LeafPositionBoundaries(i)+LeafPositionBoundaries(i+1))/2;
        if EdgeDirectionShift >= 0
        j = sum(LeafPositionBoundaries >= a & LeafPositionBoundaries <= b)-1;
        else
        j = -sum(LeafPositionBoundaries >= b & LeafPositionBoundaries <= a)+1;
        end
    if i+j <=60 && i+j >=1
        LeafJawPositions_New(i+j) = Leaf_Positions_Parallel (i);
        LeafJawPositions_New(i+j+60) = Leaf_Positions_Parallel (i+60);
    end 
end

% plot the original and new MLC locations.
LeafJawPositions_Move = transpose (LeafJawPositions_New);
figure;
y = center_value(1:60);
Leaf_A_original = LeafJawPositions(1:60);
Leaf_B_original = LeafJawPositions(61:120);
Leaf_A_final = LeafJawPositions_Move(1:60);
Leaf_B_final = LeafJawPositions_Move(61:120);
plot(Leaf_A_original,y,'og',Leaf_B_original,y,'og',Leaf_A_final,y,'*m',Leaf_B_final,y,'*m');
% axis([-200,200, -200,200]);
axis equal;

% Write the new coordinates back to a DICOM plan file
LeafJawPositions_info = dicomfind(DicomInfo,"LeafJawPositions");
LeafJawPositions_info.Value{3} = LeafJawPositions_Move;
DicomInfo_MLC_Move = dicomupdate(DicomInfo,LeafJawPositions_info);
LeafJawPositions_info2 = dicomfind(DicomInfo_MLC_Move,'LeafJawPositions');

dicomwrite([],'move.dcm',DicomInfo_MLC_Move,'CreateMode', 'copy');
% status = dicomwrite(DicomInfo_MLC_Move);
 