function [spike, v, u] = PN_neuron(dt, v, u, I)
%PN_NEURON Wrapper for the Izhikevich neurons to hold parameters.

C = 100;
a = 0.3;
b = -0.2;
c = -65;
d = 8;
k = 2;
v_r = -60;
v_t = -40;
epsilon_mean = 0;
epsilon_std = 0.05;

m = size(v,1);
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
end

