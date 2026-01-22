% function [weight_matrix_KC_EN, spike_time_PN, spike_time_KC, spike_time_EN] = GPU_network_fun(input, reward, ...
%                                         number_of_PN, number_of_KC, number_of_EN, C_I_PN_sen, C_I_PN_var, ...
%                                         connection_matrix_PN_KC, connection_matrix_KC_EN, ...
%                                         weight_matrix_PN_KC, weight_matrix_KC_EN, interval)
% Sets up a network and simulates it, run 2000ms interval.
%% Initialize parameters
% interval = 60; % [ms]
% Time step:
dt = 0.25; % [ms]
% Reward signal
BA = 0;

%% visualise the process: recordings of voltage KC, synapses PN
synapses_record = gpuArray.zeros(numPN, interval/dt);
voltage_KC = gpuArray.zeros(numKC, interval/dt);
KC_ind = 2824; % manually choose one fired KC
%%
synapses_PN_KC = gpuArray.zeros(numPN, numKC);
synapses_KC_EN = gpuArray.zeros(numKC, numEN);
spike_PN = gpuArray.zeros(numPN,1);
spike_KC = gpuArray.zeros(numKC,1);
spike_EN = gpuArray.zeros(numEN,1);
I_PN = gpuArray.zeros(numPN, 1);
I_KC = gpuArray.zeros(numKC, 1);
I_EN = gpuArray.zeros(numEN, 1);
%%
PN = [gpuArray.ones(numPN,1), gpuArray.zeros(numPN,1)]*(-60);
KC = [gpuArray.ones(numKC,1), gpuArray.zeros(numKC,1)]*(-85);
EN = [gpuArray.ones(numEN,1), gpuArray.zeros(numEN,1)]*(-60);
synaptic_tag_KC_EN = gpuArray.zeros(numKC, numEN);
concentration_BA_KC_EN = gpuArray.zeros(numKC, numEN);
t_spike_KC = gpuArray.ones(numKC,1)*(-10000);
t_spike_EN = gpuArray.ones(numEN,1)*(-10000);
% spike Recordings:
spike_time_PN = gpuArray.zeros(numPN, interval/dt);
spike_time_KC = gpuArray.zeros(numKC, interval/dt);
spike_time_EN = gpuArray.zeros(numEN, interval/dt);
pre_post_spike_occured = gpuArray.zeros(numKC, numEN);
delta_t = gpuArray.zeros(numKC, numEN);
PN_spikes = gpuArray.zeros(numPN, numKC);
KC_spikes = gpuArray.zeros(numKC, numEN);
%% Main loop - iteration over time interval/dt
for idt = 1 : interval/dt
    t = idt*dt;
%     [I_PN, BA] = train_input(t, input, reward, null_input, C_I_PN_sen, C_I_PN_var);
    if t < 10.01
        I_PN = input;
        if t == 10.0 && reward == 1
            BA = 0.5;
        end
    elseif t <= interval
        I_PN = null_input;
    end
% bsxfun speedup things!!
    % PN_KC synapses
    PN_spikes = bsxfun(@times, connection_PN_KC, spike_PN);
    synapses_PN_KC = arrayfun(@PN_KC_synapse, dt, PN_spikes, synapses_PN_KC);
    % Recordings:
    synapses_record(:,idt) = synapses_PN_KC(:,KC_ind);
    
    % KC_EN synapses
    pre_post_spike_occured = bsxfun(@max, spike_KC, spike_EN');
    KC_spikes = bsxfun(@times, connection_KC_EN, spike_KC);
    delta_t = bsxfun(@minus, t_spike_KC, t_spike_EN');
    [synapses_KC_EN, weight_matrix_KC_EN, synaptic_tag_KC_EN, concentration_BA_KC_EN] = arrayfun(@KC_EN_synapse, dt, KC_spikes, synapses_KC_EN,...
                                    weight_matrix_KC_EN, synaptic_tag_KC_EN, delta_t, pre_post_spike_occured, concentration_BA_KC_EN, BA);

    I_KC = sum(g_PN_KC*bsxfun(@times, synapses_PN_KC, (0-KC(:,1))'))';
    I_EN = sum(weight_matrix_KC_EN.*bsxfun(@times, synapses_KC_EN, (0-EN(:,1))'))';
    
    % PN neurons
    [spike_PN, PN(:,1), PN(:,2)] = arrayfun(@PN_neuron, dt, PN(:,1), PN(:,2), I_PN);
    spike_time_PN(:,idt) = spike_PN;
    % KC neurons
    [spike_KC, t_spike_KC, KC(:,1), KC(:,2)] = arrayfun(@KC_neuron, dt, t, KC(:,1), KC(:,2), I_KC, t_spike_KC);
    spike_time_KC(:,idt) = spike_KC;
    % Recordings
    voltage_KC(:,idt) = KC(:,1);
    
    % EN neurons
    [spike_EN, t_spike_EN, EN(:,1), EN(:,2)] = arrayfun(@EN_neuron, dt, t, EN(:,1), EN(:,2), I_EN, t_spike_EN);
    spike_time_EN(:,idt) = spike_EN;
end
% plot
% figure(11)
% PN_ind = find(connection_PN_KC(:,KC_ind) == 1);
% raster_PN_KC = synapses_record(PN_ind, :);
% 
% raster_PN_KC_plot = gather(raster_PN_KC);
% KC_plot = gather(voltage_KC(KC_ind,:));
% for i = 1:10 % loop through all the 10 connections PN2KC
%     subplot(12, 1, i)
%     plot(1:interval/dt, raster_PN_KC_plot(i,:), 'b');
% end
% subplot(616)
% plot(1:interval/dt, KC_plot, 'r');

% raster plot
% figure(12)
% subplot(311)
% raster_PN = repmat([1:numPN]', 1, interval/dt).*gather(spike_time_PN);
% plot(1:interval/dt, raster_PN, 'k.', 'MarkerSize', 0.1);
% subplot(312)
% raster_KC = repmat([1:numKC]', 1, interval/dt).*gather(spike_time_KC);
% plot(1:interval/dt, raster_KC, 'b.', 'MarkerSize', 0.1);
% subplot(313)
% raster_EN = repmat([1:numEN]', 1, interval/dt).*gather(spike_time_EN);
% plot(1:interval/dt, raster_EN, 'rx', 'MarkerSize', 1.0);
% axis([0 interval/dt+10, 0 3])
