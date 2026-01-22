function [S, g, c, d] = KC_EN_synapse(dt, spikes, S, g, c, delta_t, pre_post_spike_occured, d, BA )
%KC_EN_SYNAPSE Wrapper for the general synapse type to hold parameters for
%the KC-EN synapse

% tau_syn_S is the synaptic time constant.
% g is the synaptic weight or conductance (adaped in learning)
% c is the synaptic "tag"
% tau_c is a timeconstant associatec with the synaptic tag.
% delta_t = t_pre - t_post
% pre_post_spike_occured indicates if a pre- or postsynaptic spike occured.
% d is an extracellular concentration of a biogenic amine
% BA is the amount of BA released.

% Parameters
tau_c = 9; % [ms] Estimated from plot. CHANGE LATER!
tau_d = 4; % [ms] Estimated from plot. CHANGE LATER!

phi_S = 3.6; % 3 originally
tau_syn_S = 10; % [ms] 4.5 originally
g_max = 3.0; % maximal conductance
% c_min = -100.0;
[S] = synapse(dt, S, spikes, phi_S, tau_syn_S);

% [S] = arrayfun(@synapse, dt, S, spikes, phi_S, tau_syn_S);
% Change in c

dcdt = -c/tau_c;
% changes = arrayfun(@STDP, delta_t);
dcdt = dcdt + pre_post_spike_occured*STDP(delta_t);
% occured = find(pre_post_spike_occured == 1);
% dcdt(occured) = dcdt(occured) + STDP(delta_t(occured));

% Update c
% c = max(c_min, c + dcdt * dt);
c = c+dcdt*dt;
% Change in d
dddt = -d/tau_d;

% Update d
d = d + dddt * dt + BA;

% Change in g
dgdt = c * d;

% Update g
g = max(0, min(g + dgdt * dt, g_max));

end