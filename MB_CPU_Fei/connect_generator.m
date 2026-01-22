% connectivity generator

number_of_PN = 360; number_of_KC = 4000; PN2KC = 16;
connection_matrix_PN_KC = zeros(number_of_PN, number_of_KC);
baseline = 45; upperbound = 70;

for i = 1:100000
    connection_matrix_PN_KC = zeros(number_of_PN, number_of_KC);
    for ikc = 1:number_of_KC
        PN_selection = randperm(number_of_PN);
        connection_matrix_PN_KC(PN_selection(1:PN2KC),ikc) = 1;
    end
    distance = max(sum(connection_matrix_PN_KC,2)) - min(sum(connection_matrix_PN_KC,2));
    if distance <= baseline,
       save('PN_KC_16.mat', 'connection_matrix_PN_KC');
       break
    end
    
    if mod(i,1000) == 0
        i
    end
    
    if distance < upperbound
        upperbound = distance
        save(sprintf('PN_KC_F16_%d.mat', upperbound), 'connection_matrix_PN_KC');
    end
end
