clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the accuracy for interpolation
% Currently only evaluate the LJ force, equation refer to 'Efficient Calculation of Pairwise Nonbonded Forces', M. Chiu, A. Khan, M. Herbordt, FCCM2011
% Dependency: LJ_poly_interpolation_function.m (Generating the interpolation index and stored in txt file)
% Final result:
%       Fvdw_real: real result in single precision
%       Fvdw_poly: result evaluated using polynomial
%       average_error_rate: the average error rate for all the evaluated input r2
%
% By: Chen Yang
% 05/18/2018
% Boston University, CAAD Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m: interpolation order
% bin: # of bins per segment
% precision: # of datapoints for each interpolation
% min, max : range of distance

interpolation_order = 3;                % interpolation order, no larger than 3
segment_num = 12;                       % # of segment
bin_num = 256;                          % # of bins per segment
precision = 8;                          % # of datepoints when generating the polynomial index 
min_range = 0.01;                       % minimal range for the evaluation
max_range = min_range * 2^segment_num;  % maximum range for the evaluation (currently this is the cutoff radius)

cutoff = single(max_range);             % Cut-off radius
cutoff2 = cutoff * cutoff;
switchon = single(0.1);
switchon2 = single(switchon * switchon);
inv_denom = (cutoff2 - switchon2)^3;
denom = 1 / inv_denom;

r2 = single(min_range:0.01:max_range-0.01);

% initialize variables
inv_r2 = single(zeros(length(r2),1));
inv_r4 = single(zeros(length(r2),1));
inv_r6 = single(zeros(length(r2),1));
inv_r12 = single(zeros(length(r2),1));
inv_r8 = single(zeros(length(r2),1));
inv_r14 = single(zeros(length(r2),1));
s = single(zeros(length(r2),1));
ds = single(zeros(length(r2),1));
Fvdw_real = single(zeros(length(r2),1));
Fvdw_poly = single(zeros(size(r2,2),1));

% Coefficient gen (random number), independent of r
Aij = 8;
Bij = 6;

%% Evaluate the real result in single precision
for i = 1:size(r2,2)
    % smooth function
    if(r2(i) <= switchon2)
        s(i) = 1;
        ds(i) = 0;
    end
    if(r2(i) > switchon2 && r2(i) <= cutoff2)
        s(i) = (cutoff2 - r2(i)) * (cutoff2 - r2(i)) * (cutoff2 + 2*r2(i) - 3 * switchon2) * denom;
        ds(i) = 12 * (cutoff2 - r2(i)) * (switchon2 - r2(i)) * denom;
    end
    if(r2(i) > cutoff2)
        s(i) = 0;
        ds(i) = 0;
    end
    
    % calculate the real value
    inv_r2(i) = 1 / r2(i);
    inv_r4(i) = inv_r2(i) * inv_r2(i);

    inv_r6(i)  = inv_r2(i) * inv_r4(i) * s(i);
    inv_r12(i) = inv_r6(i) * inv_r6(i) * s(i);

    inv_r14(i) = inv_r12(i) * (ds(i) - 12*s(i)*inv_r2(i));
    inv_r8(i)  = inv_r6(i)  * (ds(i) -  6*s(i)*inv_r2(i));
   
    
    % calculate the VDW force
    Fvdw_real(i) = Aij * inv_r14(i) + Bij * inv_r8(i);
end

%% Generate the interpolation table (only need to run this once if the interpolation parameter remains)
%LJ_poly_interpolation_function(interpolation_order,bin_num,precision,min_range,max_range,cutoff,switchon);

%% Evaluate the interpolation result
% Load in the index data
fileID_0  = fopen('file_c0_vdw14_f_order3_bin64.txt', 'r');
fileID_1  = fopen('file_c1_vdw14_f_order3_bin64.txt', 'r');
fileID_2  = fopen('file_c2_vdw14_f_order3_bin64.txt', 'r');
fileID_3  = fopen('file_c3_vdw14_f_order3_bin64.txt', 'r');

fileID_4  = fopen('file_c0_vdw8_f_order3_bin64.txt', 'r');
fileID_5  = fopen('file_c1_vdw8_f_order3_bin64.txt', 'r');
fileID_6  = fopen('file_c2_vdw8_f_order3_bin64.txt', 'r');
fileID_7  = fopen('file_c3_vdw8_f_order3_bin64.txt', 'r');

% Fetch the index for the polynomials
read_in_c0_vdw14 = textscan(fileID_0, '%f');
read_in_c1_vdw14 = textscan(fileID_1, '%f');
read_in_c2_vdw14 = textscan(fileID_2, '%f');
read_in_c3_vdw14 = textscan(fileID_3, '%f');
read_in_c0_vdw8 = textscan(fileID_4, '%f');
read_in_c1_vdw8 = textscan(fileID_5, '%f');
read_in_c2_vdw8 = textscan(fileID_6, '%f');
read_in_c3_vdw8 = textscan(fileID_7, '%f');

% close file
fclose(fileID_0);
fclose(fileID_1);
fclose(fileID_2);
fclose(fileID_3);
fclose(fileID_4);
fclose(fileID_5);
fclose(fileID_6);
fclose(fileID_7);

% Start evaluation
for i = 1:size(r2,2)
    % Locate the segment of the current r2
    seg_ptr = 0;        % The first segment will be #0, second will be #1, etc....
    while(r2(i) >= min_range * 2^(seg_ptr+1))
        seg_ptr = seg_ptr + 1;
    end
    if(seg_ptr >= segment_num)      % if the segment pointer is larger than the maximum number of segment, then error out
        disp('Error occur: could not locate the segment for the input r2');
        exit;
    end
    
    % Locate the bin in the current segment
    segment_min = single(min_range * 2^seg_ptr);
    segment_max = segment_min * 2;
    segment_step = (segment_max - segment_min) / bin_num;
    bin_ptr = floor((r2(i) - segment_min)/segment_step) + 1;            % the bin_ptr will give which bin it locate
    
    % Calculate the index for table lookup
    lut_index = seg_ptr * bin_num + bin_ptr;
    
    % Fetch the index for the polynomials
    c0_vdw14 = single(read_in_c0_vdw14{1}(lut_index));
    c1_vdw14 = single(read_in_c1_vdw14{1}(lut_index));
    c2_vdw14 = single(read_in_c2_vdw14{1}(lut_index));
    c3_vdw14 = single(read_in_c3_vdw14{1}(lut_index));
    c0_vdw8 = single(read_in_c0_vdw8{1}(lut_index));
    c1_vdw8 = single(read_in_c1_vdw8{1}(lut_index));
    c2_vdw8 = single(read_in_c2_vdw8{1}(lut_index));
    c3_vdw8 = single(read_in_c3_vdw8{1}(lut_index));
    
    % Calculate the poly value
    vdw14 = polyval([c3_vdw14 c2_vdw14 c1_vdw14 c0_vdw14], r2(i));
    vdw8 = polyval([c3_vdw8 c2_vdw8 c1_vdw8 c0_vdw8], r2(i));
    
    % Calculate the force
    Fvdw_poly(i) = Aij * vdw14 + Bij * vdw8;
end


%% Evaluate the Error rate
diff_rate = zeros(size(r2,2),1);
for i = 1:size(r2,2)
    difference = Fvdw_poly(i) - Fvdw_real(i);
    diff_rate(i) = abs(difference / Fvdw_real(i));
end

average_error_rate = sum(diff_rate) / size(r2,2);