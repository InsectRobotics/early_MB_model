function [ change ] = STDP( delta_t )
%STDP returns the amount of change in the synaptic weight depending on the
% amount of time between pre- and postsynaptic spike.

% delta_t = t_pre - t_post

% Parameters
A_plus = 1;
A_minus = -1;
tau_plus = 3; % [ms]
tau_minus = 3; % [ms]

% change = zeros(size(delta_t));
% change = gpuArray.zeros(size(delta_t));
% gpuArray ask for element wise operation:
if delta_t < 0
    change = A_minus*exp(delta_t/tau_plus);
else
    change = A_minus*exp(-delta_t/tau_plus);
end  
% Vectorised version
% stp = find(delta_t < 0);
% std = find(delta_t >= 0);
% change = gpuArray.zeros(size(delta_t);
% change(stp) = A_minus * exp(delta_t(stp)/tau_plus);
% change(std) = A_minus * exp(-delta_t(std)/tau_minus);

end

