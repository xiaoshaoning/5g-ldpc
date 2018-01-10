% 5g ldpc encoding
% input: s, bit sequence of dimension K * 1
% output: encoded_bits, bit sequence
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

function [encoded_bits, H, Z_c] = ldpc_encode(s, base_graph_index)

K = length(s);

encoded_bits = zeros(3*K, 1);

if base_graph_index == 1
    a = 4;
    b = 22;
    c = 26;
    d = 42;
    e = 46;
    z = K/b;
    Z_c = z;
    N = 66 * Z_c;
    z = K/b;
    set_index = lifting_size_table_lookup(z);    
    load parity_check_matrices_protocol_1
    BG = parity_check_matrices_protocol_1(:, :, set_index); %#ok<NODEF>
elseif base_graph_index == 2
    a = 4;
    b = 10;
    c = 14;
    d = 38;
    e = 42;
    z = K/b;
    Z_c = z;
    N = 50 * Z_c;
    set_index = lifting_size_table_lookup(z);    
    load parity_check_matrices_protocol_2
    BG = parity_check_matrices_protocol_2(:, :, set_index); %#ok<NODEF>
else
  error('wrong base graph index in ldpc encoding.');
end

for k = (2*Z_c):(K-1)
  if s(k+1) ~= -1
    encoded_bits(k-2*Z_c+1) = s(k+1);
  else
    s(k+1) = 0;
    encoded_bits(k-2*Z_c+1) = -1;
  end    
end

% set_index = lifting_size_table_lookup(Z_c);
% BG = parity_check_matrices_protocol(:, :, set_index);

A_prime = BG(1:a, 1:b);
B_prime = BG(1:a, (b+1):c);
C_prime = BG((a+1):e, 1:b);
D_prime = BG((a+1):e, (b+1):c);

z = Z_c;

A = spalloc(a*z, b*z, nnz(A_prime + ones(size(A_prime))));

for row_index = 1:a
    for column_index = 1:b
        if A_prime(row_index, column_index) ~= -1
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(A_prime(row_index, column_index), z)+1):z, 1:mod(A_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B = spalloc(a*z, a*z, nnz(B_prime + ones(size(B_prime))));

for row_index = 1:a
    for column_index = 1:a
        if B_prime(row_index, column_index) ~= -1
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(B_prime(row_index, column_index), z)+1):z, 1:mod(B_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

C = spalloc(d*z, b*z, nnz(C_prime + ones(size(C_prime))));

for row_index = 1:d
    for column_index = 1:b
        if C_prime(row_index, column_index) ~= -1
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(C_prime(row_index, column_index), z)+1):z, 1:mod(C_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

D = spalloc(d*z, a*z, nnz(D_prime + ones(size(D_prime))));

for row_index = 1:d
    for column_index = 1:a
        if D_prime(row_index, column_index) ~= -1
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(D_prime(row_index, column_index), z)+1):z, 1:mod(D_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B_inv = spalloc(a*z, a*z, 20*z);

if (base_graph_index == 1) && (set_index ~= 7)
    B_inv(1:z, 1:z)             = speye(z);
    B_inv(1:z, 1+z:2*z)         = speye(z);
    B_inv(1:z, 1+2*z:3*z)       = speye(z);
    B_inv(1:z, 1+3*z:4*z)       = speye(z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
elseif (base_graph_index == 2) && ((set_index ~= 4) && (set_index ~= 8))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 1) && (z == 208)
    B_inv(1:z, 1:z)             = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+z:2*z)         = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+2*z:3*z)       = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+3*z:4*z)       = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1:z)         = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1:z)       = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+3*z:4*z, 1:z)       = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
elseif (base_graph_index == 1) && ((z ~= 208) && (set_index == 7))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 2) && ((set_index == 4) || (set_index == 8))    
    B_inv(1:z, 1:z)             = speye(z);
    B_inv(1:z, 1+z:2*z)         = speye(z);
    B_inv(1:z, 1+2*z:3*z)       = speye(z);
    B_inv(1:z, 1+3*z:4*z)       = speye(z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);    
end

s = s(:);

p_1 = mod(B_inv * (A * s), 2);
p_2 = mod(C * s + D * p_1, 2);

w = [p_1; p_2];

for k = K:(N+2*Z_c-1)
   encoded_bits(k-2*Z_c+1) = w(k-K+1);   
end    

H = [A, B, spalloc(a*z, d*z, 0); C, D, speye(d*z)];

% mod(H * [s; p_1; p_2]) = 0

clear A
clear B
clear C
clear D
clear B_inv
clear A_prime
clear B_prime
clear C_prime
clear D_prime

end
