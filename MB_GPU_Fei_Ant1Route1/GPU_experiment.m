% Test with the combination of 70 training images and 20 novel images
% Procedure: 
% 1.--> train with 80 images of an established route, and record activity
% 2.--> test the post-training response of the network to the 80 images
% 3.--> test the post-training response to the 5 novel images
% 4.--> record the parameter setting as well, and save the file.
% writen by Fei Peng
% 20 Sep 2013, 00:51
% reset(gpuDevice); clear all; close all; clc
tic;
% profile on
% for i_round = 1:3
%     for i_g = 1:3
%% Phase0 --> Load data and setup parameters
numPN = 360; numKC = 8000; numEN = 2; % archtecture
C_I_PN_sen = 0;     % 308 + PN-C30k1.56 works
C_I_PN_var = 3910; 
g_PN_KC = 0.25; g_KC_EN = 2.0;
numTest = 5; interval = 15; dt = 0.25;
numTrain =80;
load('PN_KC_F14_82.mat');
null_input = gpuArray.zeros(numPN,1);
connection_PN_KC = gpuArray(connection_matrix_PN_KC);
connection_KC_EN = gpuArray.ones(numKC, numEN);
% setup connection and weight matrix for KC-EN
weight_matrix_KC_EN = connection_KC_EN*g_KC_EN;
% load('./result80_win/Record80_03.mat');
% weight_matrix_KC_EN = Record.weight_matrix_KC_EN;
% Load images from raw_images folder
%**************************************************************************
%% training images
% load from ImgGrabber (workspace: struct F)
load('Raw_images_numPos_180degree80_win_evenHeading.mat'); % for training
training_inputs = zeros(numPN, numTrain);
for img = 1:numTrain,
    temp_img = Raw_images.position(img).origin;
%     temp_img = uint8((1-double(temp_img)/255)*255);
%     temp_img = imgEqual(temp_img);
    temp_img = imresize(temp_img, [10, 36]);
    temp_img = 1-double(temp_img)/255;
    temp_img = adapthisteq(temp_img);
%     temp_img = adapthisteq(temp_img);
%     temp_img = histeq(temp_img);
    temp_img = reshape(temp_img, numPN, 1);
%     subplot(5,4,img)
%     imshow(temp_img)
%     training_inputs(:,img) = reshape(temp_img, numPN, 1);
%     temp_img = (temp_img - min(temp_img))/max(temp_img)*0.8+0.1;
    training_inputs(:,img) = temp_img;
end
training_inputs = training_inputs./(repmat(sqrt(sum(training_inputs.^2)),numPN,1));
training_inputs = training_inputs*C_I_PN_var+C_I_PN_sen;
training_inputs = gpuArray(training_inputs);

%% Testing images (Control)
testing_inputs = zeros(numPN, numTest);
for img = 1:numTest,
    temp_img = Raw_images.test_position(img).origin;
%     temp_img = uint8((1-double(temp_img)/255)*255);
%     temp_img = imgEqual(temp_img);
    temp_img = imresize(temp_img, [10, 36]);
    temp_img = 1-double(temp_img)/255;
    temp_img = adapthisteq(temp_img);
    temp_img = reshape(temp_img, numPN, 1);
%     temp_img = (temp_img - min(temp_img))/max(temp_img)*0.8+0.1;
    testing_inputs(:,img) = temp_img;
%     testing_inputs(:,img) = reshape(temp_img, numPN, 1);
end
testing_inputs = testing_inputs./(repmat(sqrt(sum(testing_inputs.^2)),numPN,1));
testing_inputs = testing_inputs*C_I_PN_var+C_I_PN_sen;
testing_inputs = gpuArray(testing_inputs);
%**************************************************************************
clear Record; clear PN_activity KC_activity EN_activity

% initialize recordings in gpuArray
PN_activity = gpuArray.zeros(numTrain*2+numTest, 1);
KC_activity = gpuArray.zeros(numTrain*2+numTest, 1);
EN_activity = gpuArray.zeros(numTrain*2+numTest, 1);
KC_fired_total = gpuArray.zeros(numTrain*2+numTest, 1);
%% Phase1 --> Training phase
for n_train = 1:numTrain,
    n_train
    % reward signal
    reward = 1;
%     interval = 60;
    % call the main network file
    input = training_inputs(:,n_train);
%     weight_matrix_KC_EN = weight_matrix_KC_EN.*1.004;
    GPU_network_fun;
%     [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(training_inputs(:,n_train), reward,...
%             number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%             connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
    % main recording on each inputs
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind);
%     PN_rate = sum(spike_time_PN,2);
    KC_ind = find(sum(spike_time_KC, 2)>0);
%     KC_debug = gather(spike_time_KC);
    KC_count = length(KC_ind);
    KC_fired_total(KC_ind, n_train) = 1;
    EN_count = sum(sum(spike_time_EN));
%     EN_debug = gather(spike_time_EN);
    PN_activity(n_train) = PN_count;
    KC_activity(n_train) = KC_count;
    EN_activity(n_train) = EN_count;
end

%% Phase2 --> Post-train phase (test with the same images)
for n_train_post = 1:numTrain,
    n_train_post
    % reward signal
    reward = 0;
%     interval = 60;
    % call the main network file
    input = training_inputs(:,n_train_post);
    GPU_network_fun;
%     [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(training_inputs(:,n_train_post), reward,...
%             number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%             connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
    % main recording on each inputs
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind);
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_fired_total(KC_ind, numTrain+n_train_post) = 1;
    KC_count = length(KC_ind);
    EN_count = sum(sum(spike_time_EN));
    PN_activity(numTrain + n_train_post) = PN_count;
    KC_activity(numTrain + n_train_post) = KC_count;
    EN_activity(numTrain + n_train_post) = EN_count;
end

%% Phase3 --> Post test phase
for n_test_post = 1:numTest,
    n_test_post
    % reward signal
    reward = 0;
%     interval = 60;
    % call the main network file
    input = testing_inputs(:,n_test_post);
    GPU_network_fun;
%     [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(testing_inputs(:,n_test_post), reward,...
%             number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%             connection_matrix_PN_KC, connection_matrix_KC_EN, weight_matrix_PN_KC, weight_matrix_KC_EN, interval);
    
    % make recording on each inputs
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind);
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_fired_total(KC_ind, numTrain*2+n_test_post) = 1;
    KC_count = length(KC_ind);
    EN_count = sum(sum(spike_time_EN));
    PN_activity(numTrain * 2 + n_test_post) = PN_count;
    KC_activity(numTrain * 2 + n_test_post) = KC_count;
    EN_activity(numTrain * 2 + n_test_post) = EN_count;
    
end

%% Phase4 --> Record paras. and save file
Record.parameters.C_I_PN_sen = C_I_PN_sen; Record.parameters.C_I_PN_var = C_I_PN_var;
Record.parameters.g_PN_KC = g_PN_KC;
Record.parameters.no_PN = numPN; Record.parameters.no_KC = numKC;
Record.exp_setting = sprintf('train%d and test%d', numTrain, numTest);

% gather information from GPU
Record.PN_activity = gather(PN_activity);
Record.KC_activity = gather(KC_activity);
Record.EN_activity = gather(EN_activity);
Record.KC_fired_total = gather(KC_fired_total);
% PN_rate = gather(PN_rate);
Record.weight_matrix_KC_EN = gather(weight_matrix_KC_EN);

% close all;
figure(1)
subplot(3,1,1)
plot(1:(2*numTrain + numTest), Record.EN_activity, 'rx-');
subplot(3,1,2)
plot(1:(2*numTrain + numTest), Record.KC_activity, 'bd-');
subplot(3,1,3)
plot(1:(2*numTrain + numTest), Record.PN_activity, 'g*-');
figure(2)
hist(reshape(Record.weight_matrix_KC_EN, numKC*numEN, 1),400);
toc;
figure(3)
overlap = zeros(1, numTrain*2+numTest);
for i = 2:numTrain*2+numTest
    overlap_idx = find(KC_fired_total(:,1) == 1 & KC_fired_total(:,i) == 1);
    overlap(i) = length(overlap_idx);
    spiking_rate = sum(KC_fired_total);
end
plot(1:numTrain*2+numTest, overlap, 'k-');
hold on; plot(1:numTrain*2+numTest, spiking_rate, 'b-'); hold off;


% profile viewer
% p = profile('info');
% profsave(p, 'profile_results')