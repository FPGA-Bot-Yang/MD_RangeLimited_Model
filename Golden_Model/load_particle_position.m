% Load the particle position data from a txt file. The position data is extracted from ApoA1.pdb file.
% Unit of the position is angstrom
% The readin data is stored in a 2D array

function position_data = load_particle_position(filename)
    fp = fopen(filename);
    if fp == -1
        fprintf('failed to open %s\n',filename);
    end
    
    line_counter = 1;
    while ~feof(fp)
        tline = fgets(fp);
        line_element = textscan(tline, '%f32', 'Delimiter', '\t');         % this is cell array
        position_data(line_counter, :) = line_element{1};                  % extract the data into output array

        line_counter = line_counter + 1;
    end
end

