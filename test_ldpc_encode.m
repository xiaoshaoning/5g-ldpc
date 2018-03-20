% 5g ldpc encoding test
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT
% currently the only supported case is K = 8448

function test_ldpc_encode

K = 8448;
z = 384;

s = randi([0, 1], K, 1);

[~, ~, ~, encoded_bits_original] = ldpc_encode(s, 1);

load H

x = mod(H * encoded_bits_original, 2);

if isequal(x, zeros(46*z, 1))
    disp('LDPC encoding test passed.');
else
    disp('LDPC encoding test failed.');
end

end