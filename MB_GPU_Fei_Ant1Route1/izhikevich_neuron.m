function [spike, V, U] = izhikevich_neuron(dt, V, U, C, I, a, b, c, d, k, v_r, v_t, epsilon_mean, epsilon_std)
%IZHIKEVICH_NEURON Takes as input the current "state" of the neuron and
%outputs the estimated state after time dt.

% dt is the time step
% v is the membrane potential
% u is the recovery current
% I is the input current
% C, a, b, c, d, k, v_r, v_t, epsilon_mean, epsilon_sigma are model parameters

% m = length(V);
% C = ones(m,1)*C;
% a = ones(m,1)*a;
% b = ones(m,1)*b;
% k = ones(m,1)*k;
% v_r = ones(m,1)*v_r;
% v_t = ones(m,1)*v_t;
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
V = V + dt*((k*(V-v_r)*(V-v_t)-U+I+epsilon)/C);
% V = V + dt/2*((k*(V-v_r)*(V-v_t)-U+I+epsilon)/C);
U = U + dt*(a*(b*(V-v_r)-U));
% Reset?
if V > v_t,
V = c;
U = U + d;
spike = 1;
end


end

