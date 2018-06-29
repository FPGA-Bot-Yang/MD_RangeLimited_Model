%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Range limited force evaluation model (currently in single precision floating point)
%% Use this script to get the variable range at each step, using that range to map the data from floating point to fixed point for on-board implementation
%% 
%% By: Chen Yang
%% CAAD Lab, Boston Univerisity
%% 06/28/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Key variables:
%       position_data(particle_id, position)                    ---- the aligned position data for all the read in particles (algined to (0,0,0))
%       particle_in_cell_counter(cell_x, cell_y, cell_z)        ---- recording how many particles in each cell
%       cell_particle(cell_id, particle_id, position_data)      ---- recording the position information of each particles in each cell
%       out_range_particle_counter                              ---- the # of particles that is not in the range of the 9*9*7 cells


%% Global variables
common_path = 'F:\Research_Files\MD\Ethan_GoldenModel\Matlab_Model_Ethan\Golden_Model\';
input_position_file_name = 'input_positions_ApoA1.txt';


%% Load the data from input file
input_file_path = strcat(common_path, input_position_file_name);
fprintf('Start reading data from input file %s\n', input_file_path);
position_data = load_particle_position(input_file_path);
fprintf('Particle data loading finished!\n');
% find the min, max on each dimension
min_x  = min(position_data(:,1));
max_x  = max(position_data(:,1));
min_y  = min(position_data(:,2));
max_y  = max(position_data(:,2));
min_z  = min(position_data(:,3));
max_z  = max(position_data(:,3));
% Original range is (-56.296,56.237), (-57.123,56.259), (-40.611,40.878)
% shift all the data to positive
position_data(:,1) = position_data(:,1)-min_x;          % range: 0 ~ 112.533
position_data(:,2) = position_data(:,2)-min_y;          % range: 0 ~ 113.382
position_data(:,3) = position_data(:,3)-min_z;          % range: 0 ~ 81.489
fprintf('Particles shifted to align on (0,0,0)\n');

%% Mapping the particles to cell list
fprintf('Start mapping paricles to cells!\n');
% Bounding box of 12A, total of 9*9*7 cells, organized in a 4D array
particle_in_cell_counter = zeros(9,9,7);                        % counters tracking the # of particles in each cell
cell_particle = zeros(567,220,3);                               % 3D array holding sorted cell particles(cell_id, particle_id, position_data), cell_id = (cell_x-1)*9+(cell_y-1)*9+cell_z
out_range_particle_counter = 0;
for i = 1:92224
    % determine the cell each particle belongs to
    cell_x = ceil(position_data(i,1) / 12);
    cell_y = ceil(position_data(i,2) / 12);
    cell_z = ceil(position_data(i,3) / 12);
    % write the particle information to cell list
    if cell_x > 0 && cell_x <= 9 && cell_y > 0 && cell_y <= 9 && cell_z > 0 && cell_z <= 7
        % increment counter
        counter_temp = particle_in_cell_counter(cell_x, cell_y, cell_z) + 1;
        particle_in_cell_counter(cell_x, cell_y, cell_z) = counter_temp;
        cell_id = (cell_x-1)*9 + (cell_y-1)*9 + cell_z;
        cell_particle(cell_id,counter_temp,:) = position_data(i,:);
    else
        out_range_particle_counter = out_range_particle_counter + 1;
    end
end
fprintf('Particles mapping to cells finished! Total of %d particles falling out of the range.\n', out_range_particle_counter);
