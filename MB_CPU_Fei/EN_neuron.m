function [spike, t_spike_EN, v, u] = EN_neuron(dt, t, v, u, I, t_spike_EN)
%EN_NEURON Wrapper for the Izhikevich neurons to hold parameters.

C = 100;
a = 0.3; % original 0.3
b = -0.2; % original -0.2
c = -65; % original -65
d = 8;  % original 8
k = 2; % original 2
v_r = -60; % original -60
v_t = -40; % original -40
epsilon_mean = 0; 
epsilon_std = 0.05;

m = size(t_spike_EN,1);
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
t_spike_EN(fired) = t;

end

