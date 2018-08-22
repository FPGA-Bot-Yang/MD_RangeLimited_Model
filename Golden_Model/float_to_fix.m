%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coverting the single precision floating point to fixed point, 7.25 format, stored in uint32 integer data type, the bits are left aligned
%%
%% By: Chen Yang
%% 08/22/2018
%% Boston University, CAAD Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FINAL] = float_to_fix(in_float)

% Get the exact bits and map to an uint32 variable
X = typecast(single(in_float),'uint32');
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

end