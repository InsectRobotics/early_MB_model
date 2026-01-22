function [spike, v, u] = PN_neuron(dt, v, u, I)
%PN_NEURON Wrapper for the Izhikevich neurons to hold parameters.

C =50; %45 100
a = 0.3;
b = -0.2;
c = -65;
d =8;
k = 0.7; %2.0 %0.7
v_r = -60;
v_t = -40;
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
end
% [spike, v, u] = izhikevich_neuron(dt, v, u, C, I, a, b, c, d, k, v_r, v_t, epsilon_mean, epsilon_std);

end

