% number of trailing zeros, binary search
% according to Henry S. Warren, Jr's Hacker's Delight.
function y = ntz(x)

n = 1;

if x == 0
  y = 32;
  return
end

% if bitand(x, 65535) == 0
%   n = n + 16;
%   x = bitshift(x, -16);
% end

if bitand(x, 255) == 0
  n = n + 8;
  x = bitshift(x, -8);
end

if bitand(x, 15) == 0
  n = n + 4;
  x = bitshift(x, -4);
end

if bitand(x, 3) == 0
  n = n + 2;
  x = bitshift(x, -2);
end

y = n - bitand(x, 1);

end