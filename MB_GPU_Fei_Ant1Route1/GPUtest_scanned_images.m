% Using GPUtest_scanned_images.m to test all the scanned images I got from
% the fixed positions with a scan range of 180 degrees and 2 degree step.
tic

numScan = 91; % 91 scanned images including 1 trained
numDis = 8;
round = 2;

load('Raw_images_numPos_180degree80_win_evenHeading.mat'); % generated images
load('./result80_win/Record80_03.mat');
% load('./result80_new_PN/Record80_harder06.mat');

%% Load the paras. for network
numPN = 360; numKC = 8000; numEN = 2; % archtecture
C_I_PN_sen = 0; C_I_PN_var = 3910; %1433-1435 works
g_PN_KC = 0.25;  % g_KC_EN = 2.0;
null_input = zeros(numPN, 1);
load('PN_KC_F14_82.mat');
connection_PN_KC = gpuArray(connection_matrix_PN_KC);
connection_KC_EN = gpuArray.ones(numKC, numEN);

% setup connection and weight matrix for KC-EN
% weight_matrix_PN _KC = connection_PN_KC.*(g_PN_KC);
weight_matrix_KC_EN = gpuArray(Record.weight_matrix_KC_EN); %load the trained weights
numPos = 2; % set the number of positions want to testRaw_image
% record
KC_activity = gpuArray.zeros(numPos, numScan);
EN_activity = gpuArray.zeros(numPos, numScan);

%% Testing phase (no learning)
reward = 0; interval = 15; %[ms]

for i_pos = 1:numPos, % for each position
    i_pos
    for i_scan = 1:numScan, % for each scanned images
        i_scan
        temp_img = Raw_images.position(57+i_pos).scan_images(i_scan).raw_image;
        temp_img = imresize(temp_img, [10, 36]);
        temp_img = 1-double(temp_img)/255;
        temp_img = adapthisteq(temp_img);
        temp_img = reshape(temp_img, numPN, 1);
        test_input = temp_img;
%         temp_img = (temp_img - min(temp_img))/max(temp_img)*0.8+0.1;
        test_input = test_input./sqrt(sum(test_input.^2));
        test_input = test_input*C_I_PN_var+C_I_PN_sen;
        % call the main network file
        input = gpuArray(test_input);
        GPU_network_fun;
%         [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(test_input, reward,...
%             number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var,...
%             connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval, non_input);
%     
        % main recording on each inputs
        KC_ind = find(sum(spike_time_KC, 2)>0);
        KC_count = length(KC_ind);
        EN_count = sum(sum(spike_time_EN));
        KC_activity(i_pos, i_scan) = KC_count;
        EN_activity(i_pos, i_scan) = EN_count;
   
    end
    for i_off = 1:numDis, % for each off_loc image
        i_off
        temp_img = Raw_images.position(i_pos).off_loc(i_off).raw_image;
        temp_img = imresize(temp_img, [10, 36]);
        temp_img = 1-double(temp_img)/255;
        temp_img = adapthisteq(temp_img);
        temp_img = reshape(temp_img, numPN, 1);
        test_input = temp_img;
%         temp_img = (temp_img - min(temp_img))/max(temp_img)*0.8+0.1;
        test_input = test_input./sqrt(sum(test_input.^2));
        test_input = test_input*C_I_PN_var+C_I_PN_sen;
        input = gpuArray(test_input);
        GPU_network_fun;
%         [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(test_input, reward,...
%             number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var,...
%             connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval, non_input);
%         
        % main recording on each inputs
        KC_ind = find(sum(spike_time_KC, 2)>0);
        KC_count = length(KC_ind);
        EN_count = sum(sum(spike_time_EN));
        KC_activity(i_pos, numScan+i_off) = KC_count;
        EN_activity(i_pos, numScan+i_off) = EN_count;
    end 

end

%% Gather information and save
Record_scan_response = struct;
Record_scan_response.KC_activity = gather(KC_activity);
Record_scan_response.EN_activity = gather(EN_activity);
% save('./result80_win/Record_scan_response_03.mat', 'Record_scan_response');
% save(sprintf('./result80_new_PN/Record_scan_response%d06.mat',num_pos), 'Record_scan_response');

%% Plot all the num_pos graphs
close all;

figure(1);
for i = 1:numPos,
    no_cols = ceil(sqrt(numPos));
    no_rows = ceil(sqrt(numPos));
    subplot(no_rows, no_cols, i)
    plot(1:(numScan+numDis), Record_scan_response.EN_activity(i,:), 'rx-');
end
saveas(figure(1), sprintf('./result80_win/80posAllscanImages%d',round+1), 'jpg'); 
pause(5);

% saveas(figure(1), './result80/80pos_allscan_01','jpg');
for i = 1:numPos,
    figure(i+1)
    pause(2)
    subplot(2,1,1),
    plot(1:(numScan+numDis), Record_scan_response.EN_activity(i,:), 'rx-');
    ylabel('EN spikes')
    subplot(2,1,2),
    plot(1:(numScan+numDis), Record_scan_response.KC_activity(i,:), 'gd-');
    ylabel('KC spikes')
    hold off;
%     saveas(g,sprintf('./result80_new_PN/scan_response_P%d06',num_pos), 'jpg');
end
toc

% g = gpuDevice; reset(g);
%% do the navigation test then:
% navigate_small;