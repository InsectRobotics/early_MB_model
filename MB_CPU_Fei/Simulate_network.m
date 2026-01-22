% Sets up a network and simulates it. The first phase is the training
% phase, then comes a phase without reinforcement.
%% Initialize parameters
clear all; close all; clc
tic;
% clear all; close all; clc
Set_up_network;
% Total training time
t_end = 4000; % [ms] multiple of 1000ms
interval = 1000; 
% Time step:
dt = 0.5; % [ms]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my temp recordings for testing the net
concentration_d_1 = zeros(numKC, t_end/dt);
concentration_d_2 = zeros(numKC, t_end/dt);
synaptic_tag_c_1 = zeros(numKC, t_end/dt);
synaptic_tag_c_2 = zeros(numKC, t_end/dt);
weight_g_1 = zeros(numKC, t_end/dt);
weight_g_2 = zeros(numKC, t_end/dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reward signal
BA = 0;
% original Recordings:
% total input current per population
total_current_PN = zeros(t_end/dt,1);
total_current_KC = zeros(t_end/dt,1);
total_current_EN = zeros(t_end/dt,1);
% total_current_PCT = zeros(t_end/dt, 1);

spike_time_PN = zeros(numPN, t_end/dt);
spike_time_KC = zeros(numKC, t_end/dt);
spike_time_EN = zeros(numEN, t_end/dt);
% spike_time_PCT = zeros(number_of_PCT, t_end/dt);

% neurotransmitter per synapse population
total_neurotransmitter_PN_KC = zeros(t_end/dt,1);
total_neurotransmitter_KC_EN = zeros(t_end/dt,1);

% individual voltage per population
voltage_PN = zeros(numPN,t_end/dt);
voltage_KC = zeros(numKC,t_end/dt);
voltage_EN = zeros(numEN,t_end/dt);
% voltage_PCT = zeros(number_of_PCT, t_end/dt);

pre_post_spike_occured = zeros(numKC, numEN);
delta_t = zeros(numKC, numEN);

% Record the active rate of KC
activity_KC = zeros(t_end/dt, 1);
null_input = zeros(numPN, 1);
%% Main loop - iteration over time t_end/dt

for idt = 1 : t_end/dt
    t = idt*dt;
    if t <= interval, %interval for each training input
        [I_PN, BA] = train_input(t, inputA, 1, null_input, C_I_PN_sen, C_I_PN_var);
    elseif t <= 2*interval,
        [I_PN, BA] = train_input(t-interval, inputB, 0, null_input, C_I_PN_sen, C_I_PN_var);
    elseif t <= 3*interval,
        [I_PN, BA] = train_input(t-2*interval, inputA, 0, null_input, C_I_PN_sen, C_I_PN_var);
    elseif t <= 4*interval,
        [I_PN, BA] = train_input(t-3*interval, inputB, 0, null_input, C_I_PN_sen, C_I_PN_var);
    end
    % 2) Calculate neuro transmitter (synapses)
    
    % 2) a) PN_KC synapses   
    PN_spikes = bsxfun(@times, connection_PN_KC, spike_PN);
    synapses_PN_KC = PN_KC_synapse(dt, PN_spikes, synapses_PN_KC);
    total_neurotransmitter_PN_KC(idt) = sum(sum(synapses_PN_KC)); % Record
    
    % 2) b) KC_EN synapses
    pre_post_spike_occured = bsxfun(@max, spike_KC, spike_EN');
    KC_spikes = bsxfun(@times, connection_KC_EN, spike_KC);
    delta_t = bsxfun(@minus, t_spike_KC, t_spike_EN');
    [synapses_KC_EN, weight_matrix_KC_EN, synaptic_tag_KC_EN, concentration_BA_KC_EN] = KC_EN_synapse(dt, KC_spikes, synapses_KC_EN, ...
                                  weight_matrix_KC_EN, synaptic_tag_KC_EN, delta_t, pre_post_spike_occured, concentration_BA_KC_EN, BA);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% if the connection_KC_EN is not all to all, then the following code
%%%%% will justify the result according to the connection. But for now it
%%%%% is not necessary.
%     synapses_KC_EN = connection_matrix_KC_EN.*synapses_KC_EN;
%     weight_matrix_KC_EN = connection_matrix_KC_EN.*weight_matrix_KC_EN;
%     synaptic_tag_KC_EN = connection_matrix_KC_EN.*synaptic_tag_KC_EN;
%     concentration_BA_KC_EN = connection_matrix_KC_EN.*concentration_BA_KC_EN;
%%%%%%%%%%%%%%%%%%%%%%%%%
                              
    total_neurotransmitter_KC_EN(idt) = sum(sum(synapses_KC_EN)); % Record

    % temp record of concentration
    concentration_d_1(:,idt) = concentration_BA_KC_EN(:,1);
    concentration_d_2(:,idt) = concentration_BA_KC_EN(:,2);
    synaptic_tag_c_1(:,idt) = synaptic_tag_KC_EN(:,1);
    synaptic_tag_c_2(:,idt) = synaptic_tag_KC_EN(:,2);
    weight_g_1(:,idt) = weight_matrix_KC_EN(:,1);
    weight_g_2(:,idt) = weight_matrix_KC_EN(:,2);

    total_current_PN(idt) = sum(I_PN); % Record
    
    % 3) b) Input to KC
    I_KC = sum(g_PN_KC*bsxfun(@times, synapses_PN_KC, (0-KC(:,1))'))';
    
    % 3) c) Input to EN
    I_EN = sum(weight_matrix_KC_EN.*bsxfun(@times, synapses_KC_EN, (0-EN(:,1))'))';
    
    total_current_KC(idt) = sum(I_KC); % Record
    total_current_EN(idt) = sum(I_EN); % Record
    
    % 4) Update neurons (and record spikes)
    % 4) a) PN neurons
    [spike_PN, PN(:,1),PN(:,2)] = PN_neuron(dt, PN(:,1), PN(:,2), I_PN);
    spike_time_PN(:, idt) = spike_PN;
    voltage_PN(:,idt) = PN(:,1); % Record
    
    % 4) b) KC neurons
    [spike_KC, t_spike_KC, KC(:,1), KC(:,2)] = KC_neuron(dt, t, KC(:,1), KC(:,2), I_KC, t_spike_KC);
    spike_time_KC(:,idt) = spike_KC;
    voltage_KC(:,idt) = KC(:,1); % Record
    
    % 4) c) EN neurons
    [spike_EN, t_spike_EN, EN(:,1), EN(:,2)] = EN_neuron(dt, t, EN(:,1), EN(:,2), I_EN, t_spike_EN);
    spike_time_EN(:,idt) = spike_EN;
    voltage_EN(:,idt) = EN(:,1); % Record

end
%% plot
close all;
runtime = [runtime; toc, numKC];
figure(1);
raster_plot_4(spike_time_PN,spike_time_KC,spike_time_EN);
figure(2); plot(1:t_end/dt, concentration_d_1);
figure(3); plot(1:t_end/dt, synaptic_tag_c_1);
figure(4); plot(1:t_end/dt, weight_g_1);
figure(5);
subplot(3,2,5)
plot(total_current_PN)
title('Total current into PN');
        
hold on;
subplot(3,2,3)
plot(total_current_KC)
title('Total current into KC');
        
subplot(3,2,1)
plot(total_current_EN)
title('Total current into EN');
        
subplot(3,2,6)
plot(voltage_PN(1,:))
title('Voltage of first PN')
        
subplot(3,2,4)
plot(voltage_KC(1,:))
title('Voltage of first KC')
        
subplot(3,2,2)
plot(voltage_EN(1,:))
title('Voltage of first EN')
        
%raster_plot_2(0,0,spikes_over_time_EN);
hold off;
toc;
% saveas(figure(1), 'raster_g2410_CIN9000_A_2+.png');