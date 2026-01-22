% Ant navigation based on trained weight_matrix_KC_EN and network_fun of
% MB anti-hebbian learning form.
% 1. load data
clear all; close all; 
tic;

%% Phase 0 -- > Load data and 
load('./antview/world5000_gray.mat');
load('./antview/AntData.mat');% change here for Mike's new data
load('./result80_win/Record80_03.mat'); % Record should load to the workspace
load('Raw_images_numPos_180degree80_win_evenHeading.mat');

img_pos = Raw_images.img_pos; % load shot position from the datafile

% Details of trained route 
trained_route = img_pos;
route_length = size(trained_route,1)-1;
feeder = trained_route(1,:)';% 2D location of the feeder on XY plane
nest = trained_route(route_length,:)';% 2D location of the nest

%% Network setup
numPN = 360; numKC = 8000; numEN = 2; % archtecture
C_I_PN_sen = 0; C_I_PN_var = 3910; %1433-1435 works
g_PN_KC = 0.25; %initial_g_KC_EN = 2.0; 
interval = 15; dt = 0.25;%[ms]
reward = 0; 
load('PN_KC_F14_82.mat');
null_input = gpuArray.zeros(numPN, 1);
connection_PN_KC = gpuArray(connection_matrix_PN_KC);
% setup connection and weight matrix for KC-EN
connection_KC_EN = gpuArray.ones(numKC, numEN);
% weight_matrix_PN_KC = connection_PN_KC.*(initial_g_PN_KC);
weight_matrix_KC_EN = gpuArray(Record.weight_matrix_KC_EN);

%% TEST if the ant goes back to nest
% The infinity accounts for losting occassions
% Parameter preparation
step_size = 0.1; % [m]
step_count = 1;
infinity = 80; % full-path
% infinity = route_length; % Upper limit of steps in a navigation
%correct_threshold = 2;

scan_range = 100; % +-90 [degrees] centered in current moving direction
scan_spd = 2; % [degrees]
scan_img = scan_range/scan_spd + 1; % number of images that needed to be test
eye_height = 0.01; % [m]
resolution = 4; % [degrees/pixel]
fov = 296; % [degrees]
inputs = zeros(numPN,1);
EN_pool = zeros(1, scan_img);
%*******************************************************************
% record spikes and rotation angle at each step
% 1: absolute heading direction relative to x axis in degrees
% 2: EN spikes
% 3: if scanning at the position
step_record = zeros(3,route_length);
%*******************************************************************
% Recording positions at each step
current_position = zeros(2,route_length); 
current_position(:,1) = feeder; % Released at the feeder
% Recording the moving direction at each step
% In the first step the correct homing direction is given to the ant
% such that the ant moves one step further along the direction
moving_direction = zeros(2,route_length);

%% Phase 1 --> Scan 180 degree at step 0 to choose the right way to go
temp_scan_range = 180; temp_scan_spd = 2; initial_heading = -90; % [degree]
num_scan = temp_scan_range/temp_scan_spd + 1;
KC_activity = gpuArray.zeros(num_scan,1);
EN_activity = gpuArray.zeros(num_scan,1);
for i_scan = 1:num_scan,
    step_count
    i_scan
    temp_pos = feeder;
    temp_heading = initial_heading+temp_scan_range/2-(i_scan-1)*temp_scan_spd;
    raw_image = ImgGrabber(temp_pos(1), temp_pos(2), eye_height, temp_heading, X, Y, Z, colp, fov, resolution);
    temp_img = imresize(raw_image, [10, 36]);
    temp_img = 1-double(temp_img)/255;
    temp_img = adapthisteq(temp_img);
    temp_img = reshape(temp_img, numPN, 1);
    temp_img = temp_img./sqrt(sum(temp_img.^2));
    temp_img = temp_img*C_I_PN_var+C_I_PN_sen;
    input = gpuArray(temp_img);
    % call the main network file
    GPU_network_fun;
    % main recording on each inputs
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_count = length(KC_ind);
    EN_count = sum(sum(spike_time_EN));
    KC_activity(i_scan) = KC_count;
    EN_activity(i_scan) = EN_count;
end
% Then find the right heading to go forward:
KC_first = gather(KC_activity);
EN_first = gather(EN_activity);
[value, index] = min(EN_first);
% record heading direction and changing of moving direction (here the
% original moving direction is 0, so the first change is equal to the
% heading direction
step_record(1,1) = initial_heading+temp_scan_range/2-(index-1)*temp_scan_spd
step_record(2,1) = EN_first(index);
step_record(3,1) = 1;

%% Phase 1b --> Moving forward along current heading
moving_direction(:,step_count) = [cos(step_record(1,1)*pi/180); sin(step_record(1,1)*pi/180)];
current_position(:,step_count+1) = current_position(:,step_count) + step_size*moving_direction(:,step_count);
step_count = 2;

%% Phase 2 --> Normal procedure: (scan all the time)
for step_count = 2:infinity
    step_count
    % scan %d degrees and pre-process
    display('scanning');
    scanning = 1;   
    step_record(3,step_count) = 1;
             
    previous_heading = step_record(1,step_count-1);
    % getting all the images
    for i_scan = 1:scan_img % number of all images
        heading = previous_heading +scan_range/2 - scan_spd*(i_scan-1);
        raw_image = ImgGrabber(current_position(1,step_count),current_position(2,step_count),eye_height,heading,X,Y,Z,colp,fov,resolution);
        temp_img = imresize(raw_image, [10, 36]);
        temp_img = 1-double(temp_img)/255;
        temp_img = adapthisteq(temp_img);
        temp_img = reshape(temp_img, numPN, 1);
        temp_img = temp_img./sqrt(sum(temp_img.^2));
        temp_img = temp_img*C_I_PN_var+C_I_PN_sen;
        input = gpuArray(temp_img);
        GPU_network_fun;
        EN_count = sum(sum(spike_time_EN));
        EN_pool(i_scan) = gather(EN_count);
    end
        % Recording
        %save(sprintf('./result80/EN_pool%03d',step_count),'EN_pool');
    [step_record(2,step_count), index] = min(EN_pool);
    step_record(1,step_count) = previous_heading + scan_range/2 - scan_spd*(index-1);
    % show result
    step_record(1:2,step_count)'
    moving_direction(:,step_count) = [cos(step_record(1,step_count)*pi/180); sin(step_record(1,step_count)*pi/180)];
    current_position(:,step_count+1) = current_position(:,step_count) + step_size*moving_direction(:,step_count);
    % if reaches the nest, then break
    current_distance = norm(nest - current_position(:,step_count));
    if current_distance <= step_size
        break;
    end 
end

% plot
save('./result80_win/test_route_03','current_position');

close all

patch(X',Y',Z')
hold on
plot(trained_route(:,1),trained_route(:,2),'b','LineWidth',1.5);
plot(current_position(1,1:infinity),current_position(2,1:infinity),'r','LineWidth',0.8);

toc
        