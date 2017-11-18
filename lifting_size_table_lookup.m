% 5g ldpc parameter calculation
% input: Z
% output: set_index
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT

function set_index = lifting_size_table_lookup(Z)

for index = 1:9
  Z_shifted = bitshift(Z, -index); % floor(Z * 2^(-index));
  if Z_shifted * 2^index ~= Z
      break;
  end
end

first_element = Z * 2^(1-index);

if first_element == 1
  first_element = 2;
end

lut = [0, 1, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8];

if first_element > 15
  error('wrong input parameter Z in lift_size_table_lookup.');    
end

set_index = lut(first_element);

if set_index == 0
  error('wrong input parameter Z in lift_size_table_lookup.');    
end

end