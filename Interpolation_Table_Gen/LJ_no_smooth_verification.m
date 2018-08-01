%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using ApoA1 input postion data, to evaluate the LJ force using the generated Lookup table
% The output of this script can be used to verify the HDL simulation result
% 
% Run the following scripts before run this:
%                   LJ_no_smooth_poly_interpolation_accuracy.m          % This one generate the lookup table entries
%
% By: Chen Yang
% 07/25/2018
% Boston University, CAAD Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%% Variables
interpolation_order = 1;                % interpolation order, no larger than 3
segment_num = 14;                       % # of segment
bin_num = 256;                          % # of bins per segment
% Range starting from 2^-6 (ApoA1 min r2 is 0.015793)
min_range = 0.015625;                  % minimal range for the evaluation
% ApoA1 cutoff is 12~13 Ang, thus set the bin as 14 to cover the range
max_range = min_range * 2^segment_num;  % maximum range for the evaluation (currently this is the cutoff radius)
cutoff2 = 14*14;

filepath = '';
filename = 'input_positions_ApoA1.txt';
filename = strcat(filepath, filename);
% Position data
pos = zeros(92224,3);


%% Read in ApoA1 data
% Open File
fp = fopen(filename);
if fp == -1
    fprintf('failed to open %s\n',filename);
end
% Read in line by line
line_counter = 1;
while ~feof(fp)
    tline = fgets(fp);
    line_elements = textscan(tline,'%f');
    pos(line_counter,:) = line_elements{1}; 
    line_counter = line_counter + 1;
end
% Close File
fclose(fp);


%% Evaluate the LJ force using the generated lookup table
% Load in the index data
fileID_0  = fopen('c0_14.txt', 'r');
fileID_1  = fopen('c1_14.txt', 'r');
if interpolation_order > 1
    fileID_2  = fopen('c2_14.txt', 'r');
end
if interpolation_order > 2
    fileID_3  = fopen('c3_14.txt', 'r');
end

fileID_4  = fopen('c0_8.txt', 'r');
fileID_5  = fopen('c1_8.txt', 'r');
if interpolation_order > 1
    fileID_6  = fopen('c2_8.txt', 'r');
end
if interpolation_order > 2
    fileID_7  = fopen('c3_8.txt', 'r');
end

% Fetch the index for the polynomials
read_in_c0_vdw14 = textscan(fileID_0, '%f');
read_in_c1_vdw14 = textscan(fileID_1, '%f');
if interpolation_order > 1
    read_in_c2_vdw14 = textscan(fileID_2, '%f');
end
if interpolation_order > 2
    read_in_c3_vdw14 = textscan(fileID_3, '%f');
end
read_in_c0_vdw8 = textscan(fileID_4, '%f');
read_in_c1_vdw8 = textscan(fileID_5, '%f');
if interpolation_order > 1
    read_in_c2_vdw8 = textscan(fileID_6, '%f');
end
if interpolation_order > 2
    read_in_c3_vdw8 = textscan(fileID_7, '%f');
end

% close file
fclose(fileID_0);
fclose(fileID_1);
if interpolation_order > 1
    fclose(fileID_2);
end
if interpolation_order > 2
    fclose(fileID_3);
end
fclose(fileID_4);
fclose(fileID_5);
if interpolation_order > 1
    fclose(fileID_6);
end
if interpolation_order > 2
    fclose(fileID_7);
end

% Prepare the output file
fresult = fopen('REFERENCE_OUTPUT.txt', 'wt');
fprintf(fresult,'Total LJ\t\tX_Comp\t\tY_Comp\t\tZ_Comp\n');

pos(1,:) = [1 1 1];
pos(2,:) = [2 3 5];
pos(3,:) = [3 3 3];
pos(4,:) = [2 5 9];

% only evaluate the first reference particle
refx = single(pos(1,1));
refy = single(pos(1,2));
refz = single(pos(1,3));
for neighbor_ptr = 2:4
%for neighbor_ptr = 2:92224
    neighbor_x = single(pos(neighbor_ptr,1));
    neighbor_y = single(pos(neighbor_ptr,2));
    neighbor_z = single(pos(neighbor_ptr,3));
    
    % Calcualte r2
    dx = single(refx - neighbor_x);
    dy = single(refy - neighbor_y);
    dz = single(refz - neighbor_z);
    r2 = dx*dx + dy*dy+ dz*dz;
    % Table lookup
    % Locate the segment of the current r2
    if r2 < cutoff2
        seg_ptr = 0;        % The first segment will be #0, second will be #1, etc....
        while(r2 >= min_range * 2^(seg_ptr+1))
            seg_ptr = seg_ptr + 1;
        end
        if(seg_ptr >= segment_num)      % if the segment pointer is larger than the maximum number of segment, then error out
            disp('Error occur: could not locate the segment for the input r2');
            return;
        end
        % Locate the bin in the current segment
        segment_min = single(min_range * 2^seg_ptr);
        segment_max = single(segment_min * 2);
        segment_step = single((segment_max - segment_min) / bin_num);
        bin_ptr = floor((r2 - segment_min)/segment_step) + 1;            % the bin_ptr will give which bin it locate
        % Calculate the index for table lookup
        lut_index = seg_ptr * bin_num + bin_ptr;
        % Fetch the index for the polynomials
        c0_vdw14 = single(read_in_c0_vdw14{1}(lut_index));
        c1_vdw14 = single(read_in_c1_vdw14{1}(lut_index));
        if interpolation_order > 1
            c2_vdw14 = single(read_in_c2_vdw14{1}(lut_index));
        end
        if interpolation_order > 2
            c3_vdw14 = single(read_in_c3_vdw14{1}(lut_index));
        end
        c0_vdw8 = single(read_in_c0_vdw8{1}(lut_index));
        c1_vdw8 = single(read_in_c1_vdw8{1}(lut_index));
        if interpolation_order > 1
            c2_vdw8 = single(read_in_c2_vdw8{1}(lut_index));
        end
        if interpolation_order > 2
            c3_vdw8 = single(read_in_c3_vdw8{1}(lut_index));
        end
        % Calculate the poly value
        switch(interpolation_order)
            case 1
                vdw14 = polyval([c1_vdw14 c0_vdw14], r2);
                vdw8 = polyval([c1_vdw8 c0_vdw8], r2);
            case 2
                vdw14 = polyval([c2_vdw14 c1_vdw14 c0_vdw14], r2);
                vdw8 = polyval([c2_vdw8 c1_vdw8 c0_vdw8], r2);
            case 3
                vdw14 = polyval([c3_vdw14 c2_vdw14 c1_vdw14 c0_vdw14], r2);
                vdw8 = polyval([c3_vdw8 c2_vdw8 c1_vdw8 c0_vdw8], r2);
        end
        % Calculate the total force
        F_LJ = single(vdw14) - single(vdw8);
        F_LJ_x = single(F_LJ * dx);
        F_LJ_y = single(F_LJ * dy);
        F_LJ_z = single(F_LJ * dz);
    else
        F_LJ = 0;
        F_LJ_x = 0;
        F_LJ_y = 0;
        F_LJ_z = 0;
    end
    fprintf(fresult,'%tX\t%tX\t%tX\t%tX\n',F_LJ,F_LJ_x,F_LJ_y,F_LJ_z);
end
fclose(fresult);
