% 5g ldpc parity check matrices
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

function parity_check_matrices_calculation

parity_check_matrices_protocol_1 = zeros(46, 68, 8); 
for index = 1:8
    parity_check_matrices_protocol_1(:, :, index) = xlsread('R1-1711982_BG1.xlsx', index+1);
end

save parity_check_matrices_protocol_1 parity_check_matrices_protocol_1

parity_check_matrices_protocol_2 = zeros(42, 52, 8); 
for index = 1:8
    parity_check_matrices_protocol_2(:, :, index) = xlsread('R1-1711982_BG2.xlsx', index+1);
end

save parity_check_matrices_protocol_2 parity_check_matrices_protocol_2

end