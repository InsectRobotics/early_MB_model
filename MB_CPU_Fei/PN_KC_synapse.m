function [S] = PN_KC_synapse(dt, spikes, S)
%PN_KC_SYNAPSE Wrapper for the general synapse type to hold parameters for
%the PN-KC synapse

phi_S = 1.5;
tau_syn_S = 0.8; % [ms]

[S] = synapse(dt, S, spikes, phi_S, tau_syn_S);

end