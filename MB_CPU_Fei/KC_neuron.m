function [spike, t_spike_KC, v, u] = KC_neuron(dt, t, v, u, I, t_spike_KC)
%KC_NEURON Wrapper for the Izhikevich neurons to hold parameters.

C = 4;
a = 0.01;
b = -0.3;
c = -65;
d = 8;
k = 0.015;
v_r = -85;
v_t = -25;
epsilon_mean = 0;
epsilon_std = 0.05;

m = size(t_spike_KC,1);
spike = zeros(m,1);

% Noise term
epsilon = epsilon_mean + epsilon_std*randn(m,1);

v = v + dt*((k*(v-v_r).*(v-v_t)-u+I+epsilon)/C);
u = u + dt*(a*(b*(v-v_r)-u));
% Reset?
fired = find(v>v_t);
v(fired) = c;
u(fired) = u(fired) + d;
spike(fired) = 1;
t_spike_KC(fired) = t;
end