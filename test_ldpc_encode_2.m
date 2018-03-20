% 5g ldpc encoding test
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT
% currently the only supported case is K = 3840

function test_ldpc_encode_2

K = 3840;
z = 384;

s = randi([0, 1], K, 1);

[~, H, ~, encoded_bits_original] = ldpc_encode(s, 2);

x = mod(H * encoded_bits_original, 2);

y = zeros(42*z, 1);

if isequal(x, y)
    disp('LDPC encoding test passed.');
else
    disp('LDPC encoding test failed.');
end

end