% connectivity generator

number_of_PN = 360; number_of_KC = 8000; PN2KC = 14;
connection_matrix_PN_KC = zeros(number_of_PN, number_of_KC);
baseline = 20; upperbound = 100;

for i = 1:100000
    
    connection_matrix_PN_KC = zeros(number_of_PN, number_of_KC);
    for ikc = 1:number_of_KC
        PN_selection = randperm(number_of_PN);
        connection_matrix_PN_KC(PN_selection(1:PN2KC),ikc) = 1;
    end
    distance = max(sum(connection_matrix_PN_KC,2)) - min(sum(connection_matrix_PN_KC,2));
    if distance <= baseline,
       save('PN_KC_10.mat', 'connection_matrix_PN_KC');
    end
    
    if mod(i,1000) == 0
        i
    end
    
    if distance < upperbound
        upperbound = distance
        save(sprintf('PN_KC_F14_%d.mat', upperbound), 'connection_matrix_PN_KC');
    end
end
