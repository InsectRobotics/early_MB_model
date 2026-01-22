% Ant navigation based on trained weight_matrix_KC_EN and network_fun of
% MB anti-hebbian learning form.
% 1. load data
clear all; close all; reset(gpuDevice); 
tic;

%% Phase 0 -- > Load data and 
load('./antview/world5000_gray.mat');
load('./antview/AntData.mat');% change here for Mike's new data
load('./result20/2020record.mat'); % Record should load to the workspace

getTrainImg_Shang; % run the image generator to generate image structure and pic shot position:

% Details of trained route 
trained_route = img_pos;
route_length = size(trained_route,1);
feeder = trained_route(1,:)';% 2D location of the feeder on XY plane
nest = trained_route(route_length,:)';% 2D location of the nest

% 1) Network setup
number_of_PN = 360; number_of_KC = 4000; number_of_EN = 2; % archtecture
C_I_PN = 6180; %1433-1435 works
initial_g_PN_KC = 0.25; initial_g_KC_EN = 2.0; interval = 550; %[ms]
reward = 0; 
load('PN_KC_10op1.mat');
connection_matrix_PN_KC = gpuArray(connection_matrix_PN_KC);
% setup connection and weight matrix for KC-EN
connection_matrix_KC_EN = gpuArray.ones(number_of_KC, number_of_EN);
weight_matrix_PN_KC = connection_matrix_PN_KC.*(initial_g_PN_KC);
weight_matrix_KC_EN = connection_matrix_KC_EN.*(initial_g_KC_EN);

% 2) TEST if the ant goes back to nest
% The infinity accounts for losting occassions
% Parameter preparation
step_size = 0.1; % [m]
step_count = 2;
infinity = 20; % try with 1m
% infinity = route_length; % Upper limit of steps in a navigation
correct_threshold = 2;
scan_range = 60; % +-90 [degrees] centered in current moving direction
scan_spd = 10; % [degrees]
scan_img = scan_range/scan_spd + 1; % number of images that needed to be test
eye_height = 0.01; % [m]
resolution = 1; % [degrees/pixel]
fov = 360; % [degrees]
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
heading = -130.3464;
% moving_direction(:,1) = (nest - feeder)/norm(nest - feeder);
% current_position(:,2) = current_position(:,1) + step_size*moving_direction(:,1);
step_record(:,1) = [0, 0, 0, atan2(moving_direction(2,1),moving_direction(1,1))*180/pi];

%% Phase 1 --> Navigating

for step_count = 2:infinity
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
    temp_img = double(raw_images(1).raw_image(:,:,2))/255;
    temp_img = imresize(temp_img, [10, 36]);
    inputs = reshape(temp_img, number_of_PN, 1);
    
    % if the whole view is covered by grass,go along the original direction
    if  sum(sum(inputs)) == 0
        display('run into grass, moving forward');
        step_record(:,step_count) = [0, 0, 0, heading];
        current_position(:,step_count+1) = current_position(:,step_count)...
            + moving_direction(:,step_count)*step_size;
    else
        norm_inputs = inputs./sqrt(sum(inputs.^2));
        % norm_inputs = naive_visual_input(inputs);
        % test network with the image to get number of spikings
        display('testing')
        % call GPU_network_fun
        [weight_matrix_KC_EN, spike_time_KC, spike_time_EN] = GPU_network_fun(norm_inputs, reward, ...
                                        number_of_PN, number_of_KC, number_of_EN, C_I_PN, ...
                                        connection_matrix_PN_KC, connection_matrix_KC_EN, ...
                                        weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
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
                + moving_direction(:,step_count)*step_size;
            step_record(4,step_count) = heading;
        else
            % scan 180 degrees and pre-process
            display('unfamiliar, scanning');
            step_record(3,step_count) = 1;
            scanning = 1;            
            headings(1) = atan2(moving_direction(2,step_count),moving_direction(1,step_count))*180/pi + scan_range/2;
            % getting all the images
            for i = 1:scan_img
                heading(i) = heading(1) - scan_spd*(i-1);
                raw_images(i).raw_image = ImgGrabber(current_position(1,step_count),current_position(2,step_count),eye_height,heading(i),X,Y,Z,colp,fov,resolution);
                temp_img = double(raw_images(i).raw_image(:,:,2))/255;
                temp_img = imresize(temp_img, [10, 36]);
                inputs = reshape(temp_img, number_of_PN, 1);
                norm_inputs = inputs./sqrt(sum(inputs.^2));
                [weight_matrix_KC_EN, spike_time_KC, spike_time_EN] = GPU_network_fun(norm_inputs, reward, ...
                                        number_of_PN, number_of_KC, number_of_EN, C_I_PN, ...
                                        connection_matrix_PN_KC, connection_matrix_KC_EN, ...
                                        weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
                EN_count = sum(sum(spike_time_EN));
                EN_pool(i) = gather(EN_count);
            end
            % Recording
            save(sprintf('./result20/EN_pool%03d',step_count),'EN_pool');
            [step_record(2,step_count), index] = min(EN_pool);
            most_familiar_index = index - scan_range/(2*scan_spd) - 1;
            step_record(1,step_count) = scan_spd*most_familiar_index;
            % show result
            step_record(1:2,step_count)'
            rotation_matrix = [cos(step_record(1,step_count)/180*pi),-sin(step_record(1,step_count)/180*pi);sin(step_record(1,step_count)/180*pi),cos(step_record(1,step_count)/180*pi)];
            step_transformation = moving_direction(:,step_count)'*step_size*rotation_matrix;    
            current_position(:,step_count+1) = current_position(:,step_count) + step_transformation';
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
    if current_distance <= step_size
        break;
    end 
end

% plot
save('./result20/test_route3','current_position');

close all

patch(X',Y',Z')
hold on
plot(current_position(1,1:infinity),current_position(2,1:infinity),'r','LineWidth',1.5);
plot(trained_route(:,1),trained_route(:,2),'b','LineWidth',1.5);
toc
        