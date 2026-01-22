% Ant navigation based on trained weight_matrix_KC_EN and network_fun of
% MB anti-hebbian learning form.
% 1. load data
clear all; close all; 
tic;

%% Phase 0 -- > Load data and 
load('./antview/world5000_gray.mat');
load('./antview/AntData.mat');% change here for Mike's new data
load('./result80_new_PN/Record80.mat'); % Record should load to the workspace
load('Raw_images_numPos_180degree80.mat');

img_pos = Raw_images.img_pos; % load shot position from the datafile

% Details of trained route 
trained_route = img_pos;
small_step = 0.02; % small step
big_step = 0.1; % [m]
route_length = (size(trained_route,1)-1)/small_step*big_step*1.5; % maximum steps allowed
feeder = trained_route(1,:)';% 2D location of the feeder on XY plane
nest = trained_route(end-1,:)';% 2D location of the nest

% 1) Network setup
number_of_PN = 360; number_of_KC = 4000; number_of_EN = 2; % archtecture
C_I_PN_sen = 216.06;C_I_PN_var = 17; %1433-1435 works
initial_g_PN_KC = 0.25; initial_g_KC_EN = 2.0; interval = 550; %[ms]
reward = 0; 
load('PN_KC41.mat');
non_input = gpuArray.zeros(number_of_PN,1);
connection_matrix_PN_KC = gpuArray(connection_matrix_PN_KC);
% setup connection and weight matrix for KC-EN
connection_matrix_KC_EN = gpuArray.ones(number_of_KC, number_of_EN);
weight_matrix_PN_KC = connection_matrix_PN_KC.*(initial_g_PN_KC);
weight_matrix_KC_EN = gpuArray(Record.weight_matrix_KC_EN);

% 2) TEST if the ant goes back to nest
% The infinity accounts for losting occassions
% Parameter preparation

step_count = 1;
infinity = 80; % try with 2m
% infinity = route_length; % Upper limit of steps in a navigation
correct_threshold = 4;
hard_threshold = 50;
scan_range = 90; % +-90 [degrees] centered in current moving direction
scan_spd = 2; % [degrees]
scan_img = scan_range/scan_spd + 1; % number of images that needed to be test
eye_height = 0.01; % [m]
resolution = 4; % [degrees/pixel]
fov = 296; % [degrees]
inputs = zeros(number_of_PN,1);
EN_pool = zeros(1, scan_img);
%*******************************************************************
% record spikes and rotation angle at each step
% 1: rotation angle relative to current moving direction in degrees
% 2: EN spikes
% 3: if scanning at the position
% 4: absolute heading direction relative to x axis in degrees
step_record = zeros(4,route_length);
%*******************************************************************
% Recording positions at each step
current_position = zeros(2,route_length); 
current_position(:,1) = feeder; % Released at the feeder
% Recording the moving direction at each step
% In the first step the correct homing direction is given to the ant
% such that the ant moves one step further along the direction
moving_direction = zeros(2,route_length);
% first moving_direction comes from directly trained_route
% moving_direction(:,1) = (trained_route(2,:)' - trained_route(1,:)')/norm(trained_route(2,:)'-trained_route(1,:)');
% heading = -130.3464;
% moving_direction(:,1) = (nest - feeder)/norm(nest - feeder);
% current_position(:,2) = current_position(:,1) + step_size*moving_direction(:,1);
% step_record(:,1) = [0, 0, 0, atan2(moving_direction(2,1),moving_direction(1,1))*180/pi];

%% Phase 1 --> Scan 180 degree at step 0 to choose the right way to go
temp_scan_range = 180; temp_scan_spd = 2; initial_heading = -90; % [degree]
num_scan = temp_scan_range/temp_scan_spd + 1;
KC_activity = gpuArray.zeros(num_scan,1);
EN_activity = gpuArray.zeros(num_scan,1);
for i_scan = 1:num_scan,
    step_count
    temp_pos = feeder;
    temp_heading = initial_heading+temp_scan_range/2-(i_scan-1)*temp_scan_spd;
    raw_image = ImgGrabber(temp_pos(1), temp_pos(2), eye_height, temp_heading, X, Y, Z, colp, fov, resolution);
    temp_img = double(raw_image)/255;
    temp_img = imresize(temp_img, [10, 36]);
    temp_input = reshape(temp_img, number_of_PN, 1);
    temp_input = temp_input./sqrt(sum(temp_input.^2));
    
    % call the main network file
    input = gpuArray(temp_input);
    GPU_network_fun;
%     [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(temp_input, reward,...
%         number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var,...
%         connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
    
    % main recording on each inputs
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_count = length(KC_ind);
    EN_count = sum(sum(spike_time_EN));
    KC_activity(i_scan) = KC_count;
    EN_activity(i_scan) = EN_count;
end
% Then find the right heading to go forward:
KC_first= gather(KC_activity);
EN_first = gather(EN_activity);
save('06EN_first', 'EN_first');
[value, index] = min(EN_first);
step_record(4,1) = -90+temp_scan_range/2-(index-1)*temp_scan_spd;
step_record(1,1) = step_record(4,1)
step_record(2,1) = EN_first(index);
step_record(3,1) = 1;

%% Phase 1b --> Moving forward along current heading
moving_direction(:,step_count) = [cos(step_record(4,1)*pi/180), sin(step_record(4,1)*pi/180)];
current_position(:,step_count+1) = current_position(:,step_count) + small_step*moving_direction(:,step_count);
step_count = 2;

%% Phase 2 --> Normal procedure: move forward when confident enough with small step

for step_count = 2:route_length
    raw_images = struct;
    step_count
    scanning = 0;
    % last position and moving direction
    moving_direction(:,step_count) = (current_position(:,step_count)...
        - current_position(:,step_count-1))...
        /norm(current_position(:,step_count)...
        - current_position(:,step_count-1));
    % current heading
    heading = atan2(moving_direction(2,step_count),moving_direction(1,step_count))*180/pi;
    % get current view along current moving direction
    raw_images(1).raw_image = ImgGrabber(current_position(1,step_count),current_position(2,step_count),eye_height,heading,X,Y,Z,colp,fov,resolution);
         %************** extract feature *****************
         %inputs = extract_feature2(raw_images,step_count);
         %************** extract feature *****************
    
    % image pre-processing:
    temp_img = double(raw_images(1).raw_image)/255;
    temp_img = imresize(temp_img, [10, 36]);
    inputs = reshape(temp_img, number_of_PN, 1);
    
    % if the whole view is covered by grass,go along the original direction
    if  sum(sum(inputs)) == 0
        display('run into grass, moving forward');
        step_record(:,step_count) = [0, 0, 0, heading];
        current_position(:,step_count+1) = current_position(:,step_count)...
            + moving_direction(:,step_count)*small_step;
    else
        input = gpuArray(inputs./sqrt(sum(inputs.^2)));
        % norm_inputs = naive_visual_input(inputs);
        % test network with the image to get number of spikings
        display('testing')
        % call GPU_network_fun
        GPU_network_fun;
%         [weight_matrix_KC_EN, spike_time_KC, spike_time_EN] = GPU_network_fun(norm_inputs, reward, ...
%                                         number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%                                         connection_matrix_PN_KC, connection_matrix_KC_EN, ...
%                                         weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
        EN_count = sum(sum(spike_time_EN));
        step_record(2,step_count) = gather(EN_count);
        step_record(1,step_count) = 0; % the current forward view

%         save(sprintf('./result1/EN_count%03d',step_count),'EN_count');
%         
        % make a decision whether scan or go ahead
        if step_record(2,step_count) <= correct_threshold
            % one step further along current moving direction
            display('familiar, moving forward');
            current_position(:,step_count+1) = current_position(:,step_count)...
                + moving_direction(:,step_count)*small_step;
            step_record(4,step_count) = heading;
        else
            % scan %d degrees and pre-process
            display('unfamiliar, scanning');
            step_record(3,step_count) = 1;
            scanning = 1;            
            temp_heading = atan2(moving_direction(2,step_count),moving_direction(1,step_count))*180/pi;
            % getting all the images
            for i_scan = 1:scan_img
                heading = temp_heading + scan_range/2 - scan_spd*(i_scan-1); % scan from left to right 
                raw_images(i_scan).raw_image = ImgGrabber(current_position(1,step_count),current_position(2,step_count),eye_height,heading,X,Y,Z,colp,fov,resolution);
                temp_img = double(raw_images(i_scan).raw_image)/255;
                temp_img = imresize(temp_img, [10, 36]);
                inputs = reshape(temp_img, number_of_PN, 1);
                input = gpuArray(inputs./sqrt(sum(inputs.^2)));
                GPU_network_fun;
%                 [weight_matrix_KC_EN, spike_time_KC, spike_time_EN] = GPU_network_fun(norm_inputs, reward, ...
%                                         number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%                                         connection_matrix_PN_KC, connection_matrix_KC_EN, ...
%                                         weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
                EN_count = sum(sum(spike_time_EN));
                EN_pool(i_scan) = gather(EN_count);
            end
            % Recording
            save(sprintf('./result80_new_PN/07_EN_pool%02d',step_count),'EN_pool');
            [step_record(2,step_count), index] = min(EN_pool);
            
            if step_record(2,step_count) > hard_threshold
                current_position(:,step_count) = current_position(:,step_count-1);% MOVE BACK ONE STEP
                % and scan 180 degrees:
                temp_scan_range = 180; temp_scan_spd = 2; initial_heading = -90; % [degree]
                num_scan = temp_scan_range/temp_scan_spd + 1;
                %% Phase_backup --> Scan 180 degree at previous position to choose the right way to go
                    
                    for i_scan = 1:num_scan,
                        step_count
                        temp_pos = current_position(:,step_count);
                        temp_heading = initial_heading+temp_scan_range/2-(i_scan-1)*temp_scan_spd;
                        raw_image = ImgGrabber(temp_pos(1), temp_pos(2), eye_height, temp_heading, X, Y, Z, colp, fov, resolution);
                        temp_img = double(raw_image)/255;
                        temp_img = imresize(temp_img, [10, 36]);
                        temp_input = reshape(temp_img, number_of_PN, 1);
                        temp_input = temp_input./sqrt(sum(temp_input.^2));
    
                        % call the main network file
                        input = gpuArray(temp_input);
                        GPU_network_fun;
%     [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(temp_input, reward,...
%         number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var,...
%         connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
    
                         % main recording on each inputs
                        KC_ind = find(sum(spike_time_KC, 2)>0);
                        KC_count = length(KC_ind);
                        EN_count = sum(sum(spike_time_EN));
                        KC_activity(i_scan) = KC_count;
                        EN_activity(i_scan) = EN_count;
                    end
                    % Then find the right heading to go forward:
                        KC_back = gather(KC_activity);
                        EN_back = gather(EN_activity);
                        [value, index] = min(EN_back);
                        step_record(4,step_count) = -90+temp_scan_range/2-(index-1)*temp_scan_spd;
                        step_record(1,step_count) = step_record(4,step_count)
                        step_record(2,step_count) = EN_back(index);
                        step_record(3,step_count) = 1;
                    % moving
                        moving_direction(:,step_count) = [cos(step_record(4,1)*pi/180), sin(step_record(4,1)*pi/180)];
                        current_position(:,step_count+1) = current_position(:,step_count) + small_step*moving_direction(:,step_count);
                
            else    % if EN < hard_threshold, no need to backup, just moving toward the current heading 
                step_record(1,step_count) = scan_range/2 - scan_spd*(index-1); % if left, positive; if turn right, negative;
            % This is consistent with the atan2 way of calculating angles
%             most_familiar_index = index - scan_range/(2*scan_spd) - 1;
%             step_record(1,step_count) = scan_spd*most_familiar_index;
            % show result
                step_record(1:2,step_count)'
                
            % rotation matrix [cos(theta), -sin(theta); sin(theta),
            % cos(theta)]; counterclockwise if 0 < theta <= 90,clockwise
            % if -90 <= theta < 0.
            % Here if the angle is positive, counterclockwise. 
            % if the angle is negative(turn right), clockwise.
                rotation_matrix = [cos(step_record(1,step_count)/180*pi),-sin(step_record(1,step_count)/180*pi);sin(step_record(1,step_count)/180*pi),cos(step_record(1,step_count)/180*pi)];
                step_transformation = rotation_matrix*moving_direction(:,step_count)*small_step;
%             step_transformation = moving_direction(:,step_count)'*big_step*rotation_matrix;    
                current_position(:,step_count+1) = current_position(:,step_count) + step_transformation;
            % update the correct moving direction for the current position
                moving_direction(:,step_count) = (current_position(:,step_count+1)...
                        - current_position(:,step_count))...
                        /norm(current_position(:,step_count+1)...
                        - current_position(:,step_count));
                step_record(4,step_count) = atan2(moving_direction(2,step_count),moving_direction(1,step_count))*180/pi;
        end
    end
    % if reaches the nest, then break
    current_distance = norm(nest - current_position(:,step_count));
    if current_distance <= big_step
        break;
    end 
    end
end

% plot
save('./result80_new_PN/test_route_07','current_position');

close all

patch(X',Y',Z')
hold on
plot(current_position(1,1:step_count),current_position(2,1:step_count),'r','LineWidth',1.5);
plot(trained_route(:,1),trained_route(:,2),'b','LineWidth',1.5);
toc
        