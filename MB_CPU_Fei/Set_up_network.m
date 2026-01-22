% This script sets up a network, roughly as described in "A model of 
% non-elemental olfactory learning in Drosophila"

% Parameters:
numPN = 360;
numKC = 4000;
numEN = 2;
% number_of_PCT = 20;

C_I_PN_sen = 0; C_I_PN_var = 4000; % 250*6; 1430
g_PN_KC = 0.25; %0.175
initial_g_KC_EN = 2.0;

if ~exist('runtime')
    runtime = [];
end

%% Image load
load('Raw_images_numPos_180degree80_win.mat'); % load image file you created here
for img = 1:80
temp_img = Raw_images.position(img).origin;
temp_img = imresize(temp_img, [10, 36]);
temp_img = double(temp_img)/255;
%     figure(img),imshow(temp_img)
training_inputs(:,img) = reshape(temp_img, numPN, 1);
end
for img = 1:5
    temp_img = Raw_images.test_position(img).origin;
    temp_img = imresize(temp_img, [10 36]);
    temp_img = double(temp_img)/255;
    testing_inputs(:,img) = reshape(temp_img, numPN, 1);
end
% normalisation
training_inputs = training_inputs./repmat(sqrt(sum(training_inputs.^2)),numPN, 1);
testing_inputs = testing_inputs./repmat(sqrt(sum(testing_inputs.^2)), numPN, 1);

inputA = training_inputs(:,59);
inputB = training_inputs(:,79);

% Set up neurons [volage, current]
PN = [ones(numPN,1), zeros(numPN,1)]*(-60); % Projection neurons
KC = [ones(numKC,1), zeros(numKC,1)]*(-85); % Kenyon cells
EN = [ones(numEN,1), zeros(numEN,1)]*(-60); % Extrinsic neurons
% Set up connections
% a) PN - KC connections (PN to KC fixed whereas KC to PN varied)
load('PN_KC41.mat');
connection_PN_KC = connection_matrix_PN_KC;
% 3) b) KC - EN connections
connection_KC_EN = ones(numKC,numEN); % All KC connect to the EN.

% 4) Set up synapses

% 4) a) Initialize synaptic weight/conductance
weight_matrix_KC_EN = connection_KC_EN * initial_g_KC_EN; % FILL IN SPECIAL VALUE HERE ("low")

% 4) b) S: Amount of neuro transmitter
synapses_PN_KC = zeros(numPN,numKC);
synapses_KC_EN = zeros(numKC,numEN);
% 4) c)
synaptic_tag_KC_EN = zeros(numKC,numEN);
% 4) d)
concentration_BA_KC_EN = zeros(numKC,numEN);

% 5) Set up spike vector
spike_PN = zeros(numPN,1);
spike_KC = zeros(numKC,1);
spike_EN = zeros(numEN,1);

% 6) Set up current vector
I_PN = zeros(numPN,1);
I_KC = zeros(numKC,1);
I_EN = zeros(numEN,1);

% 7) Set up spike time vector
t_spike_KC = ones(numKC,1)*(-10000); % make sure that this does not influence the simulation.
t_spike_EN = ones(numEN,1)*(-10000); % make sure that this does not influence the simulation.