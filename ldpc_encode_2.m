% 5g ldpc encoding
% input: s, bit sequence of dimension K * 1
% output: encoded_bits, bit sequence
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

function [encoded_bits, H] = ldpc_encode_2(s)

K = length(s);

z = K/10;

set_index = lifting_size_table_lookup(z);

load parity_check_matrices_protocol_2

BG_2 = parity_check_matrices_protocol_2(:, :, set_index); %#ok<NODEF>

A_prime = BG_2(1:4, 1:10);
B_prime = BG_2(1:4, 11:14);
C_prime = BG_2(5:42, 1:10);
D_prime = BG_2(5:42, 11:14);

A = zeros(4*z, 10*z);

for row_index = 1:4
    for column_index = 1:10
        if A_prime(row_index, column_index) ~= -1
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = circshift(eye(z), -A_prime(row_index, column_index));
        else
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = zeros(z, z);
        end
    end
end

B = zeros(4*z, 4*z);

for row_index = 1:4
    for column_index = 1:4
        if B_prime(row_index, column_index) ~= -1
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = circshift(eye(z), -B_prime(row_index, column_index));
        else
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = zeros(z, z);
        end
    end
end

C = zeros(38*z, 10*z);

for row_index = 1:38
    for column_index = 1:10
        if C_prime(row_index, column_index) ~= -1
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = circshift(eye(z), -C_prime(row_index, column_index));
        else
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = zeros(z, z);
        end
    end
end

D = zeros(38*z, 4*z);

for row_index = 1:38
    for column_index = 1:4
        if D_prime(row_index, column_index) ~= -1
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = circshift(eye(z), -D_prime(row_index, column_index));
        else
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = zeros(z, z);
        end
    end
end

B_inv = zeros(4*z, 4*z);
B_inv(1:z, 1:z)             = circshift(eye(z), 1);
B_inv(1:z, 1+z:2*z)         = circshift(eye(z), 1);
B_inv(1:z, 1+2*z:3*z)       = circshift(eye(z), 1);
B_inv(1:z, 1+3*z:4*z)       = circshift(eye(z), 1);
B_inv(1+z:2*z, 1:z)         = eye(z) + circshift(eye(z), 1);
B_inv(1+z:2*z, 1+z:2*z)     = circshift(eye(z), 1);
B_inv(1+z:2*z, 1+2*z:3*z)   = circshift(eye(z), 1);
B_inv(1+z:2*z, 1+3*z:4*z)   = circshift(eye(z), 1);
B_inv(1+2*z:3*z, 1:z)       = eye(z) + circshift(eye(z), 1);
B_inv(1+2*z:3*z, 1+z:2*z)   = eye(z) + circshift(eye(z), 1);
B_inv(1+2*z:3*z, 1+2*z:3*z) = circshift(eye(z), 1);
B_inv(1+2*z:3*z, 1+3*z:4*z) = circshift(eye(z), 1);
B_inv(1+3*z:4*z, 1:z)       = circshift(eye(z), 1);
B_inv(1+3*z:4*z, 1+z:2*z)   = circshift(eye(z), 1);
B_inv(1+3*z:4*z, 1+2*z:3*z) = circshift(eye(z), 1);
B_inv(1+3*z:4*z, 1+3*z:4*z) = eye(z) + circshift(eye(z), 1);

s = s(:);

p_1 = mod(B_inv * (A * s), 2);
p_2 = mod(C * s + D * p_1, 2);

encoded_bits = [s; p_1; p_2];

H = [sparse(A), sparse(B), sparse(zeros(4*z, 38*z)); sparse(C), sparse(D), sparse(eye(38*z))];

end