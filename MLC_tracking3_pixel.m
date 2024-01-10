% Third version of MLC tracking

% V1 created a MLC squence to track tumor movement, without prediction filter.
% Conditions: Max MLC leaf speed: 25 mm/s; Control points at every 50 ms.

% V2 added gap to the MLC movement, and keeps the peripheral leaves standstill.
% V2 calculate under-exposed and over-exposed areas, geometric errors and latency

% V3 update ways to calculate under-exposed and over-exposed areas.
% This one will count the pixels.

% Jinling Zhou, 6/7/2023.

clc;
clear;
close all;

% Read leaf positions and boundaries
path = 'RP.QA303005.test.dcm';
DicomInfo = dicominfo(path);
LeafJawPositions = DicomInfo.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions;
LeafPositionBoundaries = DicomInfo.BeamSequence.Item_1.BeamLimitingDeviceSequence.Item_3.LeafPositionBoundaries;

Leaf_A_o = LeafJawPositions(1:60);
Leaf_B_o = LeafJawPositions(61:120);
N_leaf = length(Leaf_A_o);
center_value = zeros(N_leaf,1);
Leaf_width = zeros(N_leaf,1);
for i = 1: N_leaf
center_value(i) = (LeafPositionBoundaries(i)+LeafPositionBoundaries(i+1))/2;
Leaf_width(i) = LeafPositionBoundaries(i+1)-LeafPositionBoundaries(i);
end

% Targeted tumor boundaries, aligned with MLC...(for the idea that this was
% the treatment plan if the tumor was still). The targeted tumor boundaries
% will move in this exact shape with time.
Leaf_Gap = Leaf_B_o - Leaf_A_o; 
nonZeroIndex = find (Leaf_Gap);
ii = min(nonZeroIndex);
jj = max(nonZeroIndex);
Tumor_boundary_A_o = Leaf_A_o(ii:jj);
Tumor_boundary_B_o = Leaf_B_o(ii:jj);
Tumor_boundary_center_o = center_value(ii:jj);
Tumor_boundary_width = Leaf_width(ii:jj);
%Tumor area calculation
Tumor_area = sum(Tumor_boundary_width.*Leaf_Gap(ii:jj));

% Read the tumor movement data
Movement = readmatrix("sine wave_A2T3.csv");

t= Movement(:,1)*1000; % unit conversion, t in ms.
Shift_t= Movement(:,2)*10; % unit conversion, dimensions in mm. 

% moving angle with respect to edge direction
angle = pi/2;

% Only move the leaves around the tumor traces, leaving peripheral leaves motionless.
Edge_shift_max = max(Shift_t)* sin(angle);
Lower_boundary = center_value(ii)- Edge_shift_max;
Higher_boundary = center_value(jj) + Edge_shift_max;
Lower_index = ii - sum(LeafPositionBoundaries > Lower_boundary & LeafPositionBoundaries < center_value(ii));
Upper_index = jj + sum(LeafPositionBoundaries > center_value(jj) & LeafPositionBoundaries < Higher_boundary);

% Only move leaves falls within the Lower_index to the higher_index

% Target MLC position, Here, the width of leaf is taken into consideration.
% latency/speed limit of MLC is considered/added after the target is determined.
N_Tumor_move = length(t);
Target_A = zeros(N_Tumor_move,N_leaf);
Target_B = zeros(N_Tumor_move,N_leaf);
for i = 1:N_Tumor_move
    Leaf_D_Shift = Shift_t(i)* cos(angle);
    Edge_D_Shift = Shift_t(i)* sin(angle);
    % MLC target location:
    [A_f, B_f] = singlemove(Leaf_A_o,Leaf_B_o,Lower_index,Upper_index,LeafPositionBoundaries,Leaf_D_Shift,Edge_D_Shift);
    Target_A(i,:) = A_f;
    Target_B(i,:) = B_f;      
end

fileID = fopen('MLC_sequence.txt', 'a');

% MLC control point at every 50 ms, to move MLC to the most updated target location.  
i = 1;
Leaf_A_start = Leaf_A_o;
Leaf_B_start = Leaf_B_o;
k = 0;
N_control = max(t)/50 + 1;
A_control = zeros(N_control,N_leaf);
B_control = zeros(N_control,N_leaf);
center_value_t = zeros(N_control,N_leaf);
MLC_area = zeros(N_control,1);
t_space = zeros(N_control, N_leaf); 
control_trace = zeros(N_control,2);
for T_MLC = 0:50:max(t)
    k = k +1;
    if T_MLC >= t(i)
    Leaf_A_target = Target_A(i,:);
    Leaf_B_target = Target_B(i,:);
    [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,Lower_index,Upper_index,Leaf_A_target,Leaf_B_target);
    i = i+1;
    else 
    Leaf_A_target = Target_A(i-1,:);
    Leaf_B_target = Target_B(i-1,:);
    [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,Lower_index,Upper_index,Leaf_A_target,Leaf_B_target);
    end
    fprintf(fileID,'t(s)=%.3f, LeafJawPositions(mm) =', T_MLC/1000);
    fprintf(fileID, '%.2f  ', A_position,B_position);
    fprintf(fileID,'\n');
    A_control(k,:) = A_position;
    B_control(k,:) = B_position;
    center_value_t(k,:) = center_value;
    MLC_area(k) = sum((B_position - A_position).*Leaf_width);
%   for 3D plotting purposes
    for j = 1: length(center_value)
       t_space(k,j) = T_MLC/1000;
    end  
%   control trace calculation
    Trace_x_sum = 0;
    Trace_y_sum = 0;
    Trace_Count = 0;
    for index = Lower_index:Upper_index
        if B_position(index)- A_position(index) > 0.5
        Trace_Count = Trace_Count + 1;
        Trace_x_sum = center_value(index) + Trace_x_sum;
        Trace_y_sum = 0.5*(A_position(index) + B_position(index)) + Trace_y_sum;
        end
    end
    control_trace(k,1) = Trace_x_sum/Trace_Count;
    control_trace(k,2) = Trace_y_sum/Trace_Count;
%
    Leaf_A_start = A_position;
    Leaf_B_start = B_position;   
end

fclose(fileID);

% intropulated tumor location every 50ms
New_t = 0:50:max(t);
New_Shift_t = interp1(t, Shift_t, New_t); % Interpolate new tumor movement values at new time points
N_leaf_tumor = jj - ii + 1;
Tumor_boundary_A = zeros(N_control,N_leaf_tumor);
Tumor_boundary_B = zeros(N_control,N_leaf_tumor);
Tumor_boundary_center = zeros(N_control,N_leaf_tumor);
Tumor_trace = zeros(N_control,2);
% Tumor movement
for i = 1:N_control
    New_Leaf_D_Shift = New_Shift_t(i)* cos(angle);
    New_Edge_D_Shift = New_Shift_t(i)* sin(angle);
    Tumor_boundary_A(i,:) = Tumor_boundary_A_o + New_Leaf_D_Shift;
    Tumor_boundary_B(i,:) = Tumor_boundary_B_o + New_Leaf_D_Shift;
    Tumor_boundary_center(i,:) = Tumor_boundary_center_o + New_Edge_D_Shift;
    Tumor_trace(i,1) = mean(Tumor_boundary_center(i,:)); 
    Tumor_trace(i,2) = mean(0.5*(Tumor_boundary_A(i,:)+Tumor_boundary_B(i,:)));
end

% Three dimensional figure of MLC sequence
%     figure(1);
%     plot3(A_control,center_value_t,t_space,'*m',B_control,center_value_t,t_space,'*m');
%     daspect([10 10 1]);
%     hold off
 
%     figure(2);
    CNumber = 125;
% %    Relative location of tumor and MLC
%     x1 = Tumor_boundary_center(CNumber,:);
%     y1_A = Tumor_boundary_A(CNumber,:);
%     y1_B = Tumor_boundary_B(CNumber,:);
%     x2 = center_value;
%     y2_A = A_control(CNumber,:);
%     y2_B = B_control(CNumber,:);
%     plot(x1,y1_A,'og',x1,y1_B,'og',x2,y2_A,'*m',x2,y2_B,'*m');
%     axis equal;

% %    Relative location of tumor and MLC, showing the leaf width
%     figure(3);
%     x1 = [Tumor_boundary_center(CNumber,:), Tumor_boundary_center(CNumber,:)];
%     y1 = [Tumor_boundary_A(CNumber,:),Tumor_boundary_B(CNumber,:)];
%     err1 = [Tumor_boundary_width/2;Tumor_boundary_width/2];
%     x2 = [center_value; center_value];
%     y2 = [A_control(CNumber,:), B_control(CNumber,:)];
%     err2 = [Leaf_width/2;Leaf_width/2];
%     errorbar(x1,y1,err1,'horizontal','g','LineStyle','none');
%     hold on;
%     errorbar(x2,y2,err2,'horizontal','m','LineStyle','none');
%     legend('Tumor','MLC');
%     axis equal;

% area difference between MLC and tumor
MLC_lowerx = center_value - 0.5 * Leaf_width;
MLC_higherx = center_value + 0.5 * Leaf_width;
MLC_vertice_x = [MLC_lowerx; MLC_higherx];
MLC_x_A = sort(MLC_vertice_x, "ascend");
MLC_x_B = flipud(MLC_x_A);
MLC_x = [MLC_x_A; MLC_x_B];
MLC_y_A = zeros(2*N_leaf,1);
MLC_y_B = zeros(2*N_leaf,1);
MLC_y = zeros(4*N_leaf,1);
Treated_area = zeros(N_control,1);
Over_exp_area = zeros(N_control,1);
Under_exp_area = zeros(N_control,1);
% for i = 1: N_control
     i = CNumber;
    for j = 1: N_leaf
    MLC_y_A(2*j-1) = A_control(i,j);
    MLC_y_A(2*j) = A_control(i,j);
    MLC_y_B(2*j-1) = B_control(i,j);
    MLC_y_B(2*j) = B_control(i,j);
    MLC_y = [MLC_y_A; flipud(MLC_y_B)];
    end
    Tumor_center_value = Tumor_boundary_center(i,:)';
    Tumor_A = Tumor_boundary_A(i,:)';
    Tumor_B = Tumor_boundary_B(i,:)';
    [Tumor_x,Tumor_y] = vertice_sort(Tumor_center_value, Tumor_boundary_width, Tumor_A, Tumor_B);
   
    poly1 = polyshape(MLC_x,MLC_y);
    poly2 = polyshape(Tumor_x, Tumor_y);
%     plot(poly1);
%     axis equal;
%     hold on
%     plot(poly2);
%     axis equal;
    polyout = intersect(poly1,poly2);
    int_sec = polyout.Vertices;
    Xint = int_sec(:,1);
    Yint = int_sec(:,2);

    Treated_area(i) =  0.5 * sum(Xint .* circshift(Yint, 1) - Yint .* circshift(Xint, 1));

% %     [Treated_area(i)] = MonteCarlo_area(int_sec(:,1),int_sec(:,2));
% % Convert the vertices to a binary mask
% mask = poly2mask(int_sec(:,1),int_sec(:,2), imageHeight, imageWidth);
% 
% % Count the number of pixels in the masked area
% areaPixels = sum(mask(:));

    Over_exp_area(i) = MLC_area(i) - Treated_area(i);
    Under_exp_area(i) = Tumor_area - Treated_area(i);
% end  

% % Geometric error
% Geo_diff = Tumor_trace - control_trace;
% Geo_err = sqrt(sum(Geo_diff(:,1).^2 + Geo_diff(:,2).^2)/length(Geo_diff));
% Geo_err_t = zeros(N_control,1);
% % change with time
% for i = 1:N_control
%     Geo_err_t(i) = sqrt(sum(Geo_diff(1:i,1).^2 + Geo_diff(1:i,2).^2)/i); 
% end
% 
% % Latency, New_t/1000 (s) is the time stamps
% Latency = phase_delay(New_t'/1000,control_trace(:,1),Tumor_trace(:,1));

% amplitude
% Movement_control = sqrt(control_trace(:,1).^2 + control_trace(:,2).^2);
% Movement_tumor = sqrt(Tumor_trace(:,1).^2 + Tumor_trace(:,2).^2);
% figure
% plot(New_t/1000, Movement_control,"og");
% hold on
% plot(New_t/1000, Movement_tumor,"*m");


% % Tumor tracking movie
% writerObj = VideoWriter('Track_movie.avi');
% open(writerObj);
% 
% % Set up figure
% fig = figure;
% set(fig, 'Position', [100, 100, 800, 600]);
% 
% % Loop through time steps
% for m = 1:N_control
%     % Plot locations
%     x1 = Tumor_boundary_center(m,:);
%     y1_A = Tumor_boundary_A(m,:);
%     y1_B = Tumor_boundary_B(m,:);
%     x2 = center_value';
%     y2_A = A_control(m,:);
%     y2_B = B_control(m,:);
%   
%     plot([x1,x1], [y1_A,y1_B],'og',[x2,x2],[y2_A,y2_B],'*m');
%     legend({'Tumor','MLC'});
%     axis equal;
%     
%     % Set axis limits, to keep the axis still for all frames.
%     xlim([-250, 250]);
%     ylim([-70, 70]);
%     
%     % Set title and axis labels
%     title(sprintf('Time elapse = %.3f s', New_t(m)/1000));
%     xlabel('Leaf edges (mm)');
%     ylabel('MLC locations(mm)');
%     
%     % Write frame to video
%     writeVideo(writerObj, getframe(fig));
% end
% 
% % Close video writer
% close(writerObj);

% Replay movie
% implay('Track_movie.avi');

% % Single movement of bank A and bank B, no speed restrain.
% function [Leaf_A_final, Leaf_B_final] = singlemove(Leaf_A_original, Leaf_B_original, I_lower,I_upper,LeafPositionBoundaries, LeafDirectionShift, EdgeDirectionShift)
% Leaf_A_final = Leaf_A_original;
% Leaf_B_final = Leaf_B_original;
% % Move along leaf direction:
% Leaf_A_move_Parallel = Leaf_A_original + LeafDirectionShift;
% Leaf_B_move_Parallel = Leaf_B_original + LeafDirectionShift;
% Leaf_A_final(I_lower:I_upper) = Leaf_A_move_Parallel(I_lower:I_upper); 
% Leaf_B_final(I_lower:I_upper) = Leaf_B_move_Parallel(I_lower:I_upper);
% % Move along the edge direction:
% for i = I_lower: I_upper
% % Move perpendicular to leaf direciton, j is the number of edge locations
% % found within EdgeDirectionShift
%         a = (LeafPositionBoundaries(i)+LeafPositionBoundaries(i+1))/2;
%         b = a + EdgeDirectionShift;   
%      if EdgeDirectionShift >= 0
%         j = sum(LeafPositionBoundaries > a & LeafPositionBoundaries < b);
%      else
%         j = -sum(LeafPositionBoundaries > b & LeafPositionBoundaries < a);
%      end
%      if i+j <=60 && i+j >=1
%         Leaf_A_final(i+j) = Leaf_A_move_Parallel(i);
%         Leaf_B_final(i+j) = Leaf_B_move_Parallel(i);
%      end 
% end       
% end

% Move bank A and bank B with a speed restriction of 25 mm/s= 2.5 * 10^(-2) mm/ms
function [A_position, B_position] = Move_control(Leaf_A_start,Leaf_B_start,I_lower,I_upper,Leaf_A_target, Leaf_B_target)
A_position = Leaf_A_start;
B_position = Leaf_B_start;
Speed_max = 0.025; % in mm/ms
Length_max = Speed_max * 50; % 1.25 mm
Gap = Leaf_B_target - Leaf_A_target;
Gap_Initial = Leaf_B_start - Leaf_A_start;
for i = I_lower: I_upper
    % add a minimum gap of 0.5 mm to avoid leaves colliding
    if Gap(i) < 0.5 && Gap_Initial(i) >= 0.5
        Target_center = (Leaf_A_start(i) + Leaf_B_start(i)) * 0.5;
        Leaf_A_target(i) = Target_center - 0.25; 
        Leaf_B_target(i) = Target_center + 0.25;
    end

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

