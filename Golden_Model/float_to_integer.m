% Get the exact bits and map to an uint32 variable
X = typecast(single(2.368),'uint32');
% Keep the mantisa and left align
Y = bitshift(X,8);
% Set the MSB as ghost bit
Z = bitset(Y,32,1);

% Extract the expotential part
exponential_val = bitand(X, uint32(2139095040));
exponential_val = bitsrl(exponential_val, 23);
right_shift_bit = 133 - exponential_val;

% Final value: shift right the bits to align
FINAL = bitsrl(Z,right_shift_bit);