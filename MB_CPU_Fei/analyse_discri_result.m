% using data in workspace to test the overlapping of response neurons,
% the activity of KC and ENs etc.
n = t_end/2000;
A_activity = zeros(1,n); B_KC = zeros(1,n); AB_activity = zeros(1,n);
overlap_A_B = zeros(1,n); overlap_A_AB = zeros(1,n);
overlap_B_AB = zeros(1,n); overlap_A_B_AB = zeros(1,n);

fired_EN = zeros(t_end/1000, 1);

for round = 1:t_end/2000
    fired_PN_A = find(sum(spike_time_PN(:, ((2*round-2)*2000+1):(2*round-1)*2000), 2)>0);
    fired_PN_B = find(sum(spike_time_PN(:, ((2*round-1)*2000+1):2*round*2000), 2)>0);
    fired_KC_A = find(sum(spike_time_KC(:, ((2*round-2)*2000+1):(2*round-1)*2000), 2)>0);
    fired_KC_B = find(sum(spike_time_KC(:, ((2*round-1)*2000+1):2*round*2000), 2)>0);
    A_PN(round) = length(fired_PN_A);
    B_PN(round) = length(fired_PN_B);
    A_KC(round) = length(fired_KC_A);
    B_KC(round) = length(fired_KC_B); 
    
    tempA = zeros(number_of_KC, 1);tempB = zeros(number_of_KC, 1);
    
    tempA(fired_KC_A) = 1;tempB(fired_KC_B) = 1;
    
    overlap_A_B(round) = length(find(tempA == 1 & tempB == 1));
    
    fired_EN(2*round-1) = sum(sum(spike_time_EN(:, ((2*round-2)*2000+1):(2*round-1)*2000)));
    fired_EN(2*round) = sum(sum(spike_time_EN(:, ((2*round-1)*2000+1):2*round*2000)));
end
% show result
A_PN, B_PN, A_KC, B_KC, overlap_A_B
figure(10); bar(fired_EN);