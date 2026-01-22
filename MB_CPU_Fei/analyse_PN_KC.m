% analyse the KC - PN relationship

[m_KCs, n_KCs] = find(spike_time_KC == 1);
idx_KC = m_KCs(5);
idx_spike_time = n_KCs(5);
time_window = (idx_spike_time - 40) : idx_spike_time;

spike_KC = spike_time_KC(idx_KC, time_window);

idx_PN = find(connection_matrix_PN_KC(:,idx_KC) == 1);
spike_PN = spike_time_PN(idx_PN,time_window);
no_PN = length(idx_PN);
no_KC = 1;

% plot
raster_PN = repmat([1:no_PN]', 1, length(time_window)).*spike_PN;
raster_KC = repmat([(no_PN+1):(no_PN+no_KC)]', 1, length(time_window)).*spike_KC;

for i_PN = 1:no_PN,
    plot(1:length(time_window), raster_PN(i_PN,:), 'r.', 'MarkerSize', 10); hold on;
end
    plot(1:length(time_window), raster_KC, 'b.', 'MarkerSize', 10); hold off;