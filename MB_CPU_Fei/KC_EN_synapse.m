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
tau_c = 100; % [ms] Estimated from plot. CHANGE LATER!
tau_d = 15; % [ms] Estimated from plot. CHANGE LATER!

phi_S = 2.0;
tau_syn_S = 6; % [ms]
g_max = 2.0; % maximal conductance
% c_min = -0.25;
[S] = synapse(dt, S, spikes, phi_S, tau_syn_S);

% [S] = arrayfun(@synapse, dt, S, spikes, phi_S, tau_syn_S);
% Change in c

dcdt = -c/tau_c;
% changes = arrayfun(@STDP, delta_t);
dcdt = dcdt + pre_post_spike_occured.*STDP(delta_t);
% Update c
c = c + dcdt * dt;

% Change in d
dddt = -d/tau_d;

% Update d
d = d + dddt * dt + BA;

% Change in g
dgdt = c .* d;

% Update g
g = max(0, min(g + dgdt * dt, g_max));

end