% This program creates a MLC squence to track tumor movement.
% Sans prediction filter.
% Conditions: Max MLC leaf speed: 25 mm/s; Tumor moves perpendicular to the MLC leaf travel direction;
% Control points at every 50 ms.
% Jinling Zhou, 4/17/2023.
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
Leaf_width(i) = LeafPositionBoundaries(i+1)-LeafPositionBoundaries(i);
end

% Targeted tumor boundaries, aligned with MLC...(for the idea that this was
% the treatment plan if the tumor was still). The targeted tumor boundaries
% will move in this exact shape with time.
Tumor_boundary_A_o = Leaf_A_o;
Tumor_boundary_B_o = Leaf_B_o;
Tumor_boundary_center_o = center_value;
Tumor_boundary_width = Leaf_width;

% Read the tumor movement data
Movement = readmatrix("sine wave_A2T3.csv");

t= Movement(:,1)*1000; % unit conversion, t in ms.
Shift_t= Movement(:,2)*10; % unit conversion, dimensions in mm. 

% moving angle with respect to edge direction
angle = pi/2;
% Target MLC position, Here, the width of leaf is taken into consideration.
% zero latency/no speed limit on MLC
for i = 1:length(t)
    Leaf_D_Shift = Shift_t(i)* cos(angle);
    Edge_D_Shift = Shift_t(i)* sin(angle);
    % MLC target location:
    [A_f, B_f] = singlemove(Tumor_boundary_A_o,Tumor_boundary_B_o,LeafPositionBoundaries,Leaf_D_Shift,Edge_D_Shift);
    Target_A(i,:) = A_f;
    Target_B(i,:) = B_f;      
end

fileID = fopen('MLC_sequence.txt', 'a');

% MLC control point at every 50 ms, to move MLC to the most updated target location.  
i = 1;
Leaf_A_start = Leaf_A_o;
Leaf_B_start = Leaf_B_o;
k = 0;
for T_MLC = 0:50:max(t)
    k = k +1;
    if T_MLC >= t(i)
    Leaf_A_target = Target_A(i,:);
    Leaf_B_target = Target_B(i,:);
    [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,LeafPositionBoundaries,Leaf_A_target,Leaf_B_target);
    i = i+1;
    else 
    Leaf_A_target = Target_A(i-1,:);
    Leaf_B_target = Target_B(i-1,:);
    [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,LeafPositionBoundaries,Leaf_A_target,Leaf_B_target);
    end
    fprintf(fileID,'t(s)=%.3f, LeafJawPositions(mm) =', T_MLC/1000);
    fprintf(fileID, '%.2f  ', A_position,B_position);
    fprintf(fileID,'\n');
    A_control(k,:) = A_position;
    B_control(k,:) = B_position;
    center_value_t(k,:) = center_value;
    for j = 1: length(center_value)
       t_space(k,j) = T_MLC/1000;
    end  
    Leaf_A_start = A_position;
    Leaf_B_start = B_position;   
end

% intropulated tumor location at every 50ms
New_t = 0:50:max(t);
New_Shift_t = interp1(t, Shift_t, New_t); % Interpolate new tumor movement values at new time points
% Tumor movement
for i = 1:length(New_t)
    New_Leaf_D_Shift = New_Shift_t(i)* cos(angle);
    New_Edge_D_Shift = New_Shift_t(i)* sin(angle);
    Tumor_boundary_A(i,:) = Tumor_boundary_A_o + New_Leaf_D_Shift;
    Tumor_boundary_B(i,:) = Tumor_boundary_B_o + New_Leaf_D_Shift;
    Tumor_boundary_center(i,:) = Tumor_boundary_center_o + New_Edge_D_Shift;
end

% Three dimensional figure of MLC sequence
%     figure(1);
%     plot3(A_control,center_value_t,t_space,'*m',B_control,center_value_t,t_space,'*m');
%     daspect([10 10 1]);
%     hold off
 
%     figure(2);
    CNumber = 125;
   % Relative location of tumor and MLC
%     x1 = Tumor_boundary_center(CNumber,:);
%     y1_A = Tumor_boundary_A(CNumber,:);
%     y1_B = Tumor_boundary_B(CNumber,:);
%     x2 = center_value;
%     y2_A = A_control(CNumber,:);
%     y2_B = B_control(CNumber,:);
%     plot(x1,y1_A,'og',x1,y1_B,'og',x2,y2_A,'*m',x2,y2_B,'*m');
%     axis equal;

   % Relative location of tumor and MLC, showing the leaf width
    figure(3);
    x1 = [Tumor_boundary_center(CNumber,:), Tumor_boundary_center(CNumber,:)];
    y1 = [Tumor_boundary_A(CNumber,:),Tumor_boundary_B(CNumber,:)];
    err = [Leaf_width/2,Leaf_width/2];
    x2 = [center_value, center_value];
    y2 = [A_control(CNumber,:), B_control(CNumber,:)];
    errorbar(x1,y1,err,'horizontal','g','LineStyle','none');
    hold on;
    errorbar(x2,y2,err,'horizontal','m','LineStyle','none');
    legend('Tumor','MLC');
    axis equal;

fclose(fileID);

% Tumor tracking movie
writerObj = VideoWriter('Track_movie.avi');
open(writerObj);

% Set up figure
fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);

% Loop through time steps
for m = 1:length(New_t)
    % Plot locations
    x1 = Tumor_boundary_center(m,:);
    y1_A = Tumor_boundary_A(m,:);
    y1_B = Tumor_boundary_B(m,:);
    x2 = center_value;
    y2_A = A_control(m,:);
    y2_B = B_control(m,:);
  
    plot([x1,x1], [y1_A,y1_B],'og',[x2,x2],[y2_A,y2_B],'*m');
    legend({'Tumor','MLC'});
    axis equal;
    
    % Set axis limits, to keep the axis still for all frames.
    xlim([-250, 250]);
    ylim([-70, 70]);
    
    % Set title and axis labels
    title(sprintf('Time elapse = %.3f s', New_t(m)/1000));
    xlabel('Leaf edges (mm)');
    ylabel('MLC locations(mm)');
    
    % Write frame to video
    writeVideo(writerObj, getframe(fig));
end

% Close video writer
close(writerObj);

% Replay movie
% implay('Track_movie.avi');

% Single movement of bank A and bank B, no speed restrain.
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

% Move bank A and bank B with a speed restriction of 25 mm/s= 2.5 * 10^(-2) mm/ms
function [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,LeafPositionBoundaries,Leaf_A_target, Leaf_B_target)
Speed_max = 0.025; % in mm/ms
Length_max = Speed_max * 50; % 1.25 mm
for i = 1: length(LeafPositionBoundaries)-1
    if abs(Leaf_A_target(i) - Leaf_A_start(i)) > Length_max
       A_position(i) = Leaf_A_start(i) + sign(Leaf_A_target(i) - Leaf_A_start(i)) * Length_max;
    else
       A_position(i) = Leaf_A_target(i);
    end
    if abs(Leaf_B_target(i) - Leaf_B_start(i)) > Length_max
       B_position(i) = Leaf_B_start(i) + sign(Leaf_B_target(i) - Leaf_B_start(i))* Length_max;
    else
       B_position(i) = Leaf_B_target(i);
    end
end
end
