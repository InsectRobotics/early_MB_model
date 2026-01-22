function [spike] = ORN_neuron(odour, dt, afinity_vector)
%ORN_NEURON takes an odour, the time step of the simulation and the
%afinity vector as inputs and either produces a spike or not, such that
%its average spike rate will be described by the rate model given in the 
%paper "A model of non-elemental olfactory learning in Drosophila".

m = 1; % molecular Hill equivalent

odour = odour + 0.0001; % ligand concentrations should not be zero.
firing_rates = 1./(1+(afinity_vector./odour).^m);
% assuming that unit of instantaneous firing rate is: [spikes/ms]
% and that we can get the total firing rate by summing up the individual
% responses to each concentration.
% I am also wondering: In the equation above: High affinity means high
% ignorance to that ligand??

p_of_spike = sum(firing_rates, 2)*dt;
spike = zeros(size(p_of_spike));
fired = find(p_of_spike-rand(length(p_of_spike),1) > 0);
spike(fired) = 1;

end

