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
% Calculation Mode
CALCULATION_MODE = 'direct';                            % Calculation mode: 'direct' or 'table'
% Benmarck Parameters
CELL_COUNT_X = 9;
CELL_COUNT_Y = 9;
CELL_COUNT_Z = 7;
CELL_PARTICLE_MAX = 220;                                        % The maximum possible particle count in each cell
TOTAL_PARTICLE = 92224;                                         % particle count in ApoA1 benchmark
% Processing Parameters
CUTOFF_RADIUS = 12;                                             % 12 angstrom cutoff radius for ApoA1
CUTOFF_RADIUS_2 = CUTOFF_RADIUS * CUTOFF_RADIUS;                % Cutoff distance square
% MD related Parameters (source:https://github.com/pandegroup/openmm/blob/master/wrappers/python/tests/systems/test_amber_ff.xml)
EPSILON = 0.065689*0.878640;                                    % EPS(H) * EPS(O)
SIGMA = 0.247135+0.295992;                                      % SIG(H) + SIG(O)
CC = 332.0636;                                                  % Coulomb constant * Charge1 * Charge2
EWALD_COEF = 0.257952;
EWALD_COEF_2 = EWALD_COEF * EWALD_COEF;
GRAD_COEF = 0.291067;

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
particle_in_cell_counter = zeros(CELL_COUNT_X,CELL_COUNT_Y,CELL_COUNT_Z);               % counters tracking the # of particles in each cell
cell_particle = zeros(CELL_COUNT_X*CELL_COUNT_Y*CELL_COUNT_Z,CELL_PARTICLE_MAX,4);      % 3D array holding sorted cell particles(cell_id, particle_id, position_data & force value), cell_id = (cell_x-1)*9+(cell_y-1)*9+cell_z
out_range_particle_counter = 0;
for i = 1:TOTAL_PARTICLE
    % determine the cell each particle belongs to
    cell_x = ceil(position_data(i,1) / CUTOFF_RADIUS);
    cell_y = ceil(position_data(i,2) / CUTOFF_RADIUS);
    cell_z = ceil(position_data(i,3) / CUTOFF_RADIUS);
    % write the particle information to cell list
    if cell_x > 0 && cell_x <= CELL_COUNT_X && cell_y > 0 && cell_y <= CELL_COUNT_Y && cell_z > 0 && cell_z <= CELL_COUNT_Z
        % increment counter
        counter_temp = particle_in_cell_counter(cell_x, cell_y, cell_z) + 1;
        particle_in_cell_counter(cell_x, cell_y, cell_z) = counter_temp;
        cell_id = (cell_x-1)*CELL_COUNT_X + (cell_y-1)*CELL_COUNT_Y + cell_z;
        cell_particle(cell_id,counter_temp,1:3) = position_data(i,:);
    else
        out_range_particle_counter = out_range_particle_counter + 1;
    end
end
fprintf('Particles mapping to cells finished! Total of %d particles falling out of the range.\n', out_range_particle_counter);

%% Filtering
%%%%%%%%%%%%%%%%%%% In the current implementation, traverse through all the 27 neighbor cells for simple accumulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Filtering and evaluation starts!\n');
% Loopping into each cell
for home_cell_x = 1:CELL_COUNT_X
    for home_cell_y = 1:CELL_COUNT_Y
        for home_cell_z = 1:CELL_COUNT_Z
            %% Generate the 27 neighboring cell list
            % Apply the boundary conditions of neiboring cells (periodic boundary conditions)
            if home_cell_x == 1
                neighbor_cell_x_list = [CELL_COUNT_X 1 2];
            elseif home_cell_x == CELL_COUNT_X
                neighbor_cell_x_list = [CELL_COUNT_X-1 CELL_COUNT_X 1];
            else
                neighbor_cell_x_list = [home_cell_x-1 home_cell_x home_cell_x+1];
            end
            if home_cell_y == 1
                neighbor_cell_y_list = [CELL_COUNT_Y 1 2];
            elseif home_cell_y == CELL_COUNT_Y
                neighbor_cell_y_list = [CELL_COUNT_Y-1 CELL_COUNT_Y 1];
            else
                neighbor_cell_y_list = [home_cell_y-1 home_cell_y home_cell_y+1];
            end
            if home_cell_z == 1
                neighbor_cell_z_list = [CELL_COUNT_Z 1 2];
            elseif home_cell_z == CELL_COUNT_Z
                neighbor_cell_z_list = [CELL_COUNT_Z-1 CELL_COUNT_Z 1];
            else
                neighbor_cell_z_list = [home_cell_z-1 home_cell_z home_cell_z+1];
            end
            
            % Collect the home cell information
            home_cell_id = (home_cell_x-1)*CELL_COUNT_X + (home_cell_y-1)*CELL_COUNT_Y + home_cell_z;
            home_particle_counter = particle_in_cell_counter(home_cell_x,home_cell_y,home_cell_z);
            
            %% Select each particle in home cell as reference particle, traver all the particles in home cells
            for ref_particle_ptr = 1:home_particle_counter
                % Initialize the force value as 0 for each reference particle
                Force_Acc = 0;
                
                % Get reference particle coordinates
                ref_particle_pos_x = cell_particle(home_cell_id,ref_particle_ptr,1);
                ref_particle_pos_y = cell_particle(home_cell_id,ref_particle_ptr,2);
                ref_particle_pos_z = cell_particle(home_cell_id,ref_particle_ptr,3);
                
                %% Traverse 27 neighboring cells
                for neighbor_cell_x_ptr=1:3
                    for neighbor_cell_y_ptr=1:3
                        for neighbor_cell_z_ptr=1:3
                            % Get the neighbor cell id
                            neighbor_cell_id = (neighbor_cell_x_list(neighbor_cell_x_ptr)-1)*CELL_COUNT_X + (neighbor_cell_y_list(neighbor_cell_y_ptr)-1)*CELL_COUNT_Y + neighbor_cell_z_list(neighbor_cell_z_ptr);
                            neighbor_particle_counter = particle_in_cell_counter(neighbor_cell_x_list(neighbor_cell_x_ptr),neighbor_cell_y_list(neighbor_cell_y_ptr),neighbor_cell_z_list(neighbor_cell_z_ptr));
                        
                            %% Traverse all the partner particles inside the current neighbor cell
                            for neighbor_particle_ptr=1:neighbor_particle_counter
                                neighbor_particle_pos_x = cell_particle(neighbor_cell_id,neighbor_particle_ptr,1);
                                neighbor_particle_pos_y = cell_particle(neighbor_cell_id,neighbor_particle_ptr,2);
                                neighbor_particle_pos_z = cell_particle(neighbor_cell_id,neighbor_particle_ptr,3);
                                
                                % Calculate the distance square btw the reference particle and parter particle
                                dist_x_2 = (neighbor_particle_pos_x - ref_particle_pos_x)^2;
                                dist_y_2 = (neighbor_particle_pos_y - ref_particle_pos_y)^2;
                                dist_z_2 = (neighbor_particle_pos_z - ref_particle_pos_z)^2;
                                dist_2 = dist_x_2 + dist_y_2 + dist_z_2;
                                
                                % Filtering, and direct force evaluation
                                A = 12 * EPSILON *(SIGMA^12);
                                B = 6 * EPSILON *(SIGMA^6);
                                if dist_2 <= CUTOFF_RADIUS_2 && dist_2 > 0
                                    inv_dist_2 = 1 / dist_2;
                                    inv_dist_4 = inv_dist_2 ^ 2;
                                    inv_dist_8 = inv_dist_4 ^ 2;
                                    inv_dist_14 = inv_dist_4 * inv_dist_8 * inv_dist_2;
                                    %% Force evaluation (??????????????????? Is the formular here right ??????????????????)
                                    % Direct computation
                                    if strcmp(CALCULATION_MODE, 'direct')
                                        Coulomb_Force = CC * erfc(EWALD_COEF*sqrt(dist_2)) * (1/sqrt(dist_2));
                                        LJ_Force = A * inv_dist_14 - B * inv_dist_8;
                                    % Table look-up
                                    elseif strcmp(CALCULATION_MODE, 'table')
                                     
                                    else
                                        fprintf('Please select a valid force evaluation mode: direct or table...\n');
                                        exit;
                                    end
                                    
                                    % Accumulate the force
                                    Force_Acc = Force_Acc + Coulomb_Force + LJ_Force;
                                end
                            end
                        end
                    end
                end
                % Assign the total force to cell_particle
                cell_particle(home_cell_id,ref_particle_ptr,4) = Force_Acc;
            end
            
        end
    end
end
