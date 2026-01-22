function [spike, t_spike_KC, v, u] = KC_neuron(dt, t, v, u, I, t_spike_KC)
%KC_NEURON Wrapper for the Izhikevich neurons to hold parameters.

C = 1.8;  %1.9
a = 0.01;
b = -0.3;
c = -65;
d = 8;
k = 0.05; % 0.015
v_r = -85;
v_t = -25;
epsilon_mean = 0;
epsilon_std = 0.05;

spike = 0;

% Noise term
epsilon = epsilon_mean + epsilon_std * randn;

% Change terms
% dvdt = (k*(V - v_r)*(V - v_t) - U + I + epsilon)/C;
% dudt = a*(b*(V-v_r) - U);
% 
% % After change
% V = V + dvdt * dt;
% U = U + dudt * dt;
v = v + dt*((k*(v-v_r)*(v-v_t)-u+I+epsilon)/C);
% V = V + dt/2*((k*(V-v_r)*(V-v_t)-U+I+epsilon)/C);
u = u + dt*(a*(b*(v-v_r)-u));
% Reset?
if v > v_t,
    v = c;
    u = u + d;
    spike = 1;
    t_spike_KC = t;
end
% [spike, v, u] = izhikevich_neuron(dt, v, u, C, I, a, b, c, d, k, v_r, v_t, epsilon_mean, epsilon_std);

end