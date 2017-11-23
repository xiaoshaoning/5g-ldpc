% 5g ldpc encoding
% input: s, bit sequence of dimension K * 1
% output: encoded_bits, bit sequence
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

function encoded_bits = ldpc_encode(s)

K = length(s);

z = K/22;

set_index = lifting_size_table_lookup(z);

load parity_check_matrices_protocol_1

BG_1 = parity_check_matrices_protocol_1(:, :, set_index); %#ok<NODEF>

A_prime = BG_1(1:4, 1:22);
B_prime = BG_1(1:4, 23:26);
C_prime = BG_1(5:46, 1:22);
D_prime = BG_1(5:46, 23:26);

A = spalloc(4*z, 22*z, nnz(A_prime + ones(size(A_prime))));

for row_index = 1:4
    for column_index = 1:22
        if A_prime(row_index, column_index) ~= -1
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [A_prime(row_index, column_index)+1:z, 1:A_prime(row_index, column_index)], ones(1, z), z, z);
        else
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B = spalloc(4*z, 4*z, nnz(B_prime + ones(size(B_prime))));

for row_index = 1:4
    for column_index = 1:4
        if B_prime(row_index, column_index) ~= -1
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [B_prime(row_index, column_index)+1:z, 1:B_prime(row_index, column_index)], ones(1, z), z, z);
        else
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

C = spalloc(42*z, 22*z, nnz(C_prime + ones(size(C_prime))));

for row_index = 1:42
    for column_index = 1:22
        if C_prime(row_index, column_index) ~= -1
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [C_prime(row_index, column_index)+1:z, 1:C_prime(row_index, column_index)], ones(1, z), z, z);
        else
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

D = spalloc(42*z, 4*z, nnz(D_prime + ones(size(D_prime))));

for row_index = 1:42
    for column_index = 1:4
        if D_prime(row_index, column_index) ~= -1
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [D_prime(row_index, column_index)+1:z, 1:D_prime(row_index, column_index)], ones(1, z), z, z);
        else
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B_inv = spalloc(4*z, 4*z, 20*z);
B_inv(1:z, 1:z) = sparse(eye(z));
B_inv(1:z, 1+z:2*z) = sparse(eye(z));
B_inv(1:z, 1+2*z:3*z) = sparse(eye(z));
B_inv(1:z, 1+3*z:4*z) = sparse(eye(z));
B_inv(1+z:2*z, 1:z) = sparse(eye(z) + circshift(eye(z), -1));
B_inv(1+z:2*z, 1+z:2*z) = sparse(circshift(eye(z), -1));
B_inv(1+z:2*z, 1+2*z:3*z) = sparse(circshift(eye(z), -1));
B_inv(1+z:2*z, 1+3*z:4*z) = sparse(circshift(eye(z), -1));
B_inv(1+2*z:3*z, 1:z) = sparse(circshift(eye(z), -1));
B_inv(1+2*z:3*z, 1+z:2*z) = sparse(circshift(eye(z), -1));
B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(eye(z) + circshift(eye(z), -1));
B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(eye(z) + circshift(eye(z), -1));
B_inv(1+3*z:4*z, 1:z) = sparse(circshift(eye(z), -1));
B_inv(1+3*z:4*z, 1+z:2*z) = sparse(circshift(eye(z), -1));
B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(circshift(eye(z), -1));
B_inv(1+3*z:4*z, 1+3*z:4*z) = sparse(eye(z) + circshift(eye(z), -1));

s = s(:);

p_1 = mod(B_inv * (A * s), 2);
p_2 = mod(C * s + D * p_1, 2);

encoded_bits = [s; p_1; p_2];

end