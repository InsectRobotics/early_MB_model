function [S] = PN_KC_synapse(dt, spikes, S)
% Parameters

phi_S = 4.7; %3.3
tau_syn_S = 0.26; % [ms]
[S] = synapse(dt, S, spikes, phi_S, tau_syn_S);

end