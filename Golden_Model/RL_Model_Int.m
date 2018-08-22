%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Range limited force evaluation model
%% (The implementation is complicated and currently paused, please check the fixed point version in file 'RL_Model_Fixed.m')
%% All data in custom fixed point type
%% Dependency:
%%                      float_to_fix
%%
%%
%% By: Chen Yang
%% CAAD Lab, Boston Univerisity
%% 08/22/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Key variables:
%       position_data(particle_id, position)                    ---- the aligned position data for all the read in particles (algined to (0,0,0))
%       particle_in_cell_counter(cell_x, cell_y, cell_z)        ---- recording how many particles in each cell
%       cell_particle(cell_id, particle_id, position_data)      ---- recording the position information of each particles in each cell
%       out_range_particle_counter                              ---- the # of particles that is not in the range of the 9*9*7 cells

clear all;

%% Global variables
ENABLE_VERIFICATION = false;
RANGE_PROFILING = true;
% Input Path
common_path = 'F:\Research_Files\MD\Ethan_GoldenModel\Matlab_Model_Ethan\Golden_Model\';
input_position_file_name = 'input_positions_ApoA1.txt';
% Calculation Mode
CALCULATION_MODE = 'direct';                                    % Select Calculation mode between 'direct' or 'table'
% Select Force Model
FORCE_MODEL = 'OpenMM';                                         % Select force model between 'OpenMM' or 'CAAD'
% Benmarck Parameters
CELL_COUNT_X = 16;
CELL_COUNT_Y = 16;
CELL_COUNT_Z = 11;
CELL_PARTICLE_MAX = 220;                                        % The maximum possible particle count in each cell
TOTAL_PARTICLE = 92224;                                         % particle count in ApoA1 benchmark
% Processing Parameters
BOND_DISTANCE = single(0.96);                                   % Bonded pairs distance for H2O molecule
BOND_DISTANCE_2 = BOND_DISTANCE ^ 2;
CUTOFF_RADIUS = single(8);                                     % 12 angstrom cutoff radius for ApoA1
CUTOFF_RADIUS_2 = CUTOFF_RADIUS * CUTOFF_RADIUS;                % Cutoff distance square
INV_CUTOFF_RADIUS = 1 / CUTOFF_RADIUS;
INV_CUTOFF_RADIUS_3 = 1 / (CUTOFF_RADIUS_2 * CUTOFF_RADIUS);    % Cutoff distance cube inverse
SWITCH_DIST = 10;                                               % Switch distance for LJ evaluation
% MD related Parameters (source:https://github.com/pandegroup/openmm/blob/master/wrappers/python/tests/systems/test_amber_ff.xml)
EPSILON = single(sqrt(0.065689*0.878640));                      % sqrt[EPS(H) * EPS(O)]
SIGMA = single((0.247135+0.295992)/2);                          % [SIG(H) + SIG(O)] / 2
CC = single(332.0636);                                          % Coulomb constant * Charge1 * Charge2
EWALD_COEF = single(0.257952);
EWALD_COEF_2 = EWALD_COEF * EWALD_COEF;
GRAD_COEF = 0.291067;
ONE_4PI_EPS0 = 138.935456;                                      % Coulomb constant, unit kJ/nm
SOLVENT_DIELECTRIC = 80;                                        % Dielectric constant of the solvent, in water, value is 80
Q1 = 1;                                                         % Charge 1 (serve as placeholder)
Q2 = 2;                                                         % Charge 2 (serve as placeholder)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Preprocessing Input data (in single precision)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Load the data from input file
% input_file_path = strcat(common_path, input_position_file_name);
% fprintf('*** Start reading data from input file %s ***\n', input_file_path);
% raw_position_data = single(load_raw_particle_position(input_file_path));
% fprintf('Particle data loading finished!\n');
% % find the min, max on each dimension
% min_x  = min(raw_position_data(:,1));
% max_x  = max(raw_position_data(:,1));
% min_y  = min(raw_position_data(:,2));
% max_y  = max(raw_position_data(:,2));
% min_z  = min(raw_position_data(:,3));
% max_z  = max(raw_position_data(:,3));
% % Original range is (-56.296,56.237), (-57.123,56.259), (-40.611,40.878)
% % shift all the data to positive
% raw_position_data(:,1) = raw_position_data(:,1)-min_x;          % range: 0 ~ 112.533
% raw_position_data(:,2) = raw_position_data(:,2)-min_y;          % range: 0 ~ 113.382
% raw_position_data(:,3) = raw_position_data(:,3)-min_z;          % range: 0 ~ 81.489
% fprintf('All particles shifted to align on (0,0,0)\n');
% 
% %% Changing the position data from single to custom fixed point in uint32 type (7.25 format, 32-bit)
% fprintf('*** Changing the format from single precision to custom fixed point!\n');
% for i = 1:TOTAL_PARTICLE
%     position_data(i,1) = float_to_fix(raw_position_data(i,1), 7);
%     position_data(i,2) = float_to_fix(raw_position_data(i,2), 7);
%     position_data(i,3) = float_to_fix(raw_position_data(i,3), 7);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the pre-processed data in custom fixed point (7.25 format)
input_file_path = strcat(common_path, input_position_file_name);
fprintf('*** Start reading data from input file %s ***\n', input_file_path);
% Load the data into workspace
load(input_file_path);
fprintf('Particle data loading finished!\n');


%% Mapping the particles to cell list
fprintf('*** Start mapping paricles to cells! ***\n');
% % Bounding box of 8A, total of 16*16*11 cells, organized in a 4D array
% % The 4 MSB determines the cell
% particle_in_cell_counter = zeros(CELL_COUNT_X,CELL_COUNT_Y,CELL_COUNT_Z);               % counters tracking the # of particles in each cell
% cell_particle = uint32(zeros(CELL_COUNT_X*CELL_COUNT_Y*CELL_COUNT_Z,CELL_PARTICLE_MAX,8));      % 3D array holding sorted cell particles(cell_id, particle_id, particle_info), cell_id = (cell_x-1)*9*7+(cell_y-1)*7+cell_z
%                                                                                         % Patticle info: 1~3:position, 4~6:force component in each direction, 7: energy, 8:# of partner particles
% out_range_particle_counter = 0;
% for i = 1:TOTAL_PARTICLE
%     % Position data shift right 32-4=28bits, then it's the cell # belongs to
%     cell_x = bitsrl(position_data(i,1),28) + 1;
%     cell_y = bitsrl(position_data(i,2),28) + 1;
%     cell_z = bitsrl(position_data(i,3),28) + 1;
%     % write the particle information to cell list
%     if cell_x > 0 && cell_x <= CELL_COUNT_X && cell_y > 0 && cell_y <= CELL_COUNT_Y && cell_z > 0 && cell_z <= CELL_COUNT_Z
%         % increment counter
%         counter_temp = particle_in_cell_counter(cell_x, cell_y, cell_z) + 1;
%         particle_in_cell_counter(cell_x, cell_y, cell_z) = counter_temp;
%         cell_id = (cell_x-1)*CELL_COUNT_Y*CELL_COUNT_Z + (cell_y-1)*CELL_COUNT_Z + cell_z;
%         cell_particle(cell_id,counter_temp,1) = position_data(i,1);
%         cell_particle(cell_id,counter_temp,2) = position_data(i,2);
%         cell_particle(cell_id,counter_temp,3) = position_data(i,3);
%     else
%         out_range_particle_counter = out_range_particle_counter + 1;
%         fprintf('Out of range partcile is (%f,%f,%f)\n', position_data(i,1:3));
%     end
% end
% fprintf('Particles mapping to cells finished! Total of %d particles falling out of the range.\n', out_range_particle_counter);

% Load the data into workspace
load('cell_particle.mat');
load('particle_in_cell_counter.mat');
fprintf('Particles mapping to cells finished!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering & Force Evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% In the current implementation, traverse through all the 27 neighbor cells for simple accumulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('*** Filtering and evaluation starts! ***\n');
fprintf('Force Model is %s, calculation mode is: %s!\n', FORCE_MODEL, CALCULATION_MODE);
switch FORCE_MODEL
    case 'OpenMM'
        fprintf('Starting evaluation using the OpenMM force model! (http://docs.openmm.org/6.2.0/userguide/theory.html)\n');
    case 'CAAD'
        fprintf('Starting evaluation using the CAAD Lab force model! (https://www.bu.edu/caadlab/hprcta10_tc.pdf)\n');
    otherwise
        fprintf('Invalid force model!\n');
end
force_write_back_counter = 0;       % Profiling variable to record how many neighbor particles has been calculated
% Loopping into each home cell
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
            home_cell_id = (home_cell_x-1)*CELL_COUNT_Y*CELL_COUNT_Z + (home_cell_y-1)*CELL_COUNT_Z + home_cell_z;
            home_particle_counter = particle_in_cell_counter(home_cell_x,home_cell_y,home_cell_z);
            
            %% Select each particle in home cell as reference particle, traverse all the particles in home cells
            for ref_particle_ptr = 1:home_particle_counter
                % Initialize the force and energy value as 0 for each reference particle
                Force_Acc_x = single(0);
                Force_Acc_y = single(0);
                Force_Acc_z = single(0);
                Energy_Acc = single(0);
                % Initialize the counter for recording how many neighboring particles within the cutoff radius
                particles_within_cutoff = 0;
                
                % Get reference particle coordinates
                ref_particle_pos_x = cell_particle(home_cell_id,ref_particle_ptr,1);
                ref_particle_pos_y = cell_particle(home_cell_id,ref_particle_ptr,2);
                ref_particle_pos_z = cell_particle(home_cell_id,ref_particle_ptr,3);
                
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % DEBUGGING VARIABLE
%                 if home_cell_x == 3 && home_cell_y == 4 && home_cell_z == 3 && ref_particle_ptr == 1
%                     counter_particle_within_cutoff = zeros(CELL_COUNT_X, CELL_COUNT_Y, CELL_COUNT_Z);
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %% Traverse 27 neighboring cells
                for neighbor_cell_x_ptr=1:3
                    for neighbor_cell_y_ptr=1:3
                        for neighbor_cell_z_ptr=1:3
                            % Get the neighbor cell id
                            neighbor_cell_x = neighbor_cell_x_list(neighbor_cell_x_ptr);
                            neighbor_cell_y = neighbor_cell_y_list(neighbor_cell_y_ptr);
                            neighbor_cell_z = neighbor_cell_z_list(neighbor_cell_z_ptr);
                            neighbor_cell_id = (neighbor_cell_x-1)*CELL_COUNT_Y*CELL_COUNT_Z + (neighbor_cell_y-1)*CELL_COUNT_Z + neighbor_cell_z;
                            neighbor_particle_counter = particle_in_cell_counter(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z);
                        
                            %% Traverse all the partner particles inside the current neighbor cell
                            for neighbor_particle_ptr=1:neighbor_particle_counter
                                neighbor_particle_pos_x = cell_particle(neighbor_cell_id,neighbor_particle_ptr,1);
                                neighbor_particle_pos_y = cell_particle(neighbor_cell_id,neighbor_particle_ptr,2);
                                neighbor_particle_pos_z = cell_particle(neighbor_cell_id,neighbor_particle_ptr,3);
                                
                                % Apply boundary condition to the partner particles
                                if home_cell_x == CELL_COUNT_X && neighbor_cell_x == 1
                                    neighbor_particle_pos_x = single(neighbor_particle_pos_x + CELL_COUNT_X * CUTOFF_RADIUS);
                                elseif home_cell_x == 1 && neighbor_cell_x == CELL_COUNT_X
                                    neighbor_particle_pos_x = single(neighbor_particle_pos_x - CELL_COUNT_X * CUTOFF_RADIUS);
                                end
                                if home_cell_y == CELL_COUNT_Y && neighbor_cell_y == 1
                                    neighbor_particle_pos_y = single(neighbor_particle_pos_y + CELL_COUNT_Y * CUTOFF_RADIUS);
                                elseif home_cell_y == 1 && neighbor_cell_y == CELL_COUNT_Y
                                    neighbor_particle_pos_y = single(neighbor_particle_pos_y - CELL_COUNT_Y * CUTOFF_RADIUS);
                                end
                                if home_cell_z == CELL_COUNT_Z && neighbor_cell_z == 1
                                    neighbor_particle_pos_z = single(neighbor_particle_pos_z + CELL_COUNT_Z * CUTOFF_RADIUS);
                                elseif home_cell_z == 1 && neighbor_cell_z == CELL_COUNT_Z
                                    neighbor_particle_pos_z = single(neighbor_particle_pos_z - CELL_COUNT_Z * CUTOFF_RADIUS);
                                end
                                
                                % Calculate the distance square btw the reference particle and parter particle
                                dist_x = single(neighbor_particle_pos_x - ref_particle_pos_x);
                                dist_y = single(neighbor_particle_pos_y - ref_particle_pos_y);
                                dist_z = single(neighbor_particle_pos_z - ref_particle_pos_z);
                                dist_x_2 = dist_x^2;
                                dist_y_2 = dist_y^2;
                                dist_z_2 = dist_z^2;
                                dist_2 = dist_x_2 + dist_y_2 + dist_z_2;
                                dist = sqrt(dist_2);
                                inv_dist = 1 / dist;
                                
                                % RANGE PROFILING
                                if RANGE_PROFILING
                                    min_dx = min([min_dx dist_x dist_y dist_z]);
                                    max_dx = max([max_dx dist_x dist_y dist_z]);
                                    min_dx_2 = min([min_dx_2 dist_x_2 dist_y_2 dist_z_2]);
                                    max_dx_2 = max([max_dx_2 dist_x_2 dist_y_2 dist_z_2]);
                                    min_r2 = min([min_r2 dist_2]);
                                    max_r2 = max([max_r2 dist_2]);
                                end
                                
                                %% Filtering, and direct force evaluation
                                % PROFILING: Count the # of bonded particle pairs
                                if dist_2 <= BOND_DISTANCE_2
                                    Bonded_Particle_Pairs = Bonded_Particle_Pairs + 1;
                                end
                                if dist_2 <= CUTOFF_RADIUS_2 && dist_2 > BOND_DISTANCE_2
                                    % increment the counter
                                    particles_within_cutoff = particles_within_cutoff + 1;

%                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     % DEBUGGING FRAGMENT
%                                     if home_cell_x == 3 && home_cell_y == 4 && home_cell_z == 3 && ref_particle_ptr == 20
%                                         counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z) = counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z) + 1;
%                                         if neighbor_cell_x == 3 && neighbor_cell_y == 4 && neighbor_cell_z == 2
%                                             temp_cell_particles_within_cutoff(counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z),1:3) = [neighbor_particle_pos_x, neighbor_particle_pos_y, neighbor_particle_pos_z];
%                                             %fprintf("%f,%f,%f\n", neighbor_particle_pos_x, neighbor_particle_pos_y, neighbor_particle_pos_z);
%                                         end
%                                     end
%                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                    % Select the force model
                                    switch FORCE_MODEL
                                        % Follow the force model in OpenMM:
                                        % https://github.com/pandegroup/openmm/blob/master/tests/TestNonbondedForce.h
                                        % Applied cutoff and switch function, but no PME
                                        case 'OpenMM'
                                            % Calculate the switch function value (http://docs.openmm.org/6.2.0/userguide/theory.html)
                                            x = (dist - SWITCH_DIST) / (CUTOFF_RADIUS - SWITCH_DIST);
                                            if dist >= SWITCH_DIST && dist <= CUTOFF_RADIUS
                                                Switch_Val = 1 - 10*x^3 + 15*x^4 - 6*x^5;
                                                Switch_Deri = (-30*x^2 + 60*x^3 - 30*x^4) / (CUTOFF_RADIUS - SWITCH_DIST);
                                            elseif dist >= 0 && dist < SWITCH_DIST
                                                Switch_Val = 1;
                                                Switch_Deri = 0;
                                            else
                                                Switch_Val = 0;
                                                Switch_Deri = 0;
                                            end
                                            %% Force evaluation
                                            sigma_6 = SIGMA^6;
                                            sigma_12 = sigma_6^2;
                                            inv_dist_2 = 1 / dist_2;
                                            inv_dist_3 = inv_dist_2 / dist;
                                            inv_dist_6 = inv_dist_2 ^ 3;
                                            inv_dist_12 = inv_dist_6 ^ 2;
                                            inv_dist_8 = inv_dist_6 * inv_dist_2;
                                            inv_dist_14 = inv_dist_6 * inv_dist_8;
                                            % Coulomb interaction with cutoff using reaction field approximation
                                            krf = INV_CUTOFF_RADIUS_3 * (SOLVENT_DIELECTRIC - 1) / (2 * SOLVENT_DIELECTRIC + 1);
                                            crf = INV_CUTOFF_RADIUS * (3 * SOLVENT_DIELECTRIC) / (2 * SOLVENT_DIELECTRIC + 1);
                                            % Direct computation
                                            if strcmp(CALCULATION_MODE, 'direct')
                                                %% Energy
                                                % LJ Potential
                                                LJ_Energy = single(Switch_Val * 4 * (EPSILON*sigma_12*inv_dist_12 - EPSILON*sigma_6*inv_dist_6));
                                                % Coulomb Potential
                                                chargeProd = ONE_4PI_EPS0 * Q1 * Q2;
                                                Coulomb_Energy = single(chargeProd * (inv_dist + krf * dist_2 - crf));
                                                % Total energy
                                                Total_Energy = single(LJ_Energy + Coulomb_Energy);
                                                %% Force (The force calculate here is F/r, for easy calculation of force component on each axis)
                                                % LJ force over R
                                                LJ_Force_over_R = Switch_Val*4*EPSILON*(12*sigma_12*inv_dist_14 - 6*sigma_6*inv_dist_8);
                                                % Apply switch condition for LJ force
                                                LJ_Force_over_R = LJ_Force_over_R - Switch_Deri*4*EPSILON*(sigma_12*inv_dist_12 - sigma_6*inv_dist_6)*inv_dist;
                                                % Coulomb Force over R
                                                Coulomb_Force_over_R = chargeProd * (inv_dist_3 - 2*krf);
                                                % Total force over R
                                                Total_Force_over_R = LJ_Force_over_R + Coulomb_Force_over_R;
                                            % Table look-up
                                            elseif strcmp(CALCULATION_MODE, 'table')
                                                fprintf('Table lookup version for OpenMM force model is under construction......\n');
                                                return;
                                            else
                                                fprintf('Please select a valid force evaluation mode: direct or table...\n');
                                                return;
                                            end

                                            % Accumulate the force & energy
                                            Energy_Acc = Energy_Acc + Total_Energy;
                                            % Total force component in each direction
                                            Total_Force_x = Total_Force_over_R * dist_x;
                                            Total_Force_y = Total_Force_over_R * dist_y;
                                            Total_Force_z = Total_Force_over_R * dist_z;
                                            % Accumulate force in each direction
                                            Force_Acc_x = Force_Acc_x + Total_Force_x;
                                            Force_Acc_y = Force_Acc_y + Total_Force_y;
                                            Force_Acc_z = Force_Acc_z + Total_Force_z;

                                        % Follow the force model from CAAD lab publications:
                                        % https://ieeexplore.ieee.org/document/5771251/
                                        % https://ieeexplore.ieee.org/document/5670800/
                                        % Applied cutoff and switch function, and PME
                                        case 'CAAD'
                                            A = 48 * EPSILON *(SIGMA^12);
                                            B = 24 * EPSILON *(SIGMA^6);
                                            % Force calculate (??????????????????? Is the formular here right ??????????????????)
                                            inv_dist_2 = 1 / dist_2;
                                            inv_dist_4 = inv_dist_2 ^ 2;
                                            inv_dist_8 = inv_dist_4 ^ 2;
                                            inv_dist_14 = inv_dist_4 * inv_dist_8 * inv_dist_2;

                                            %% Force evaluation
                                            % Direct computation
                                            if strcmp(CALCULATION_MODE, 'direct')
                                                Coulomb_Force_over_R = CC * erfc(EWALD_COEF*sqrt(dist_2)) * (1/sqrt(dist_2));
                                                LJ_Force_over_R = A * inv_dist_14 - B * inv_dist_8;
                                                Total_Force_over_R = Coulomb_Force_over_R + LJ_Force_over_R;
                                            % Table look-up
                                            elseif strcmp(CALCULATION_MODE, 'table')
                                                fprintf('Table lookup version for CAAD force model is under construction......\n');
                                                return;
                                            else
                                                fprintf('Please select a valid force evaluation mode: direct or table...\n');
                                                return;
                                            end

                                            % Accumulate the force
                                            Force_Acc_x = Force_Acc_x + Total_Force_over_R * dist_x;
                                            Force_Acc_y = Force_Acc_y + Total_Force_over_R * dist_y;
                                            Force_Acc_z = Force_Acc_z + Total_Force_over_R * dist_z;

                                        otherwise
                                            fprintf('Please select a valid force model between OpenMM or CAAD!\n');
                                            return;
                                    end
                                    
                                    % RANGE PROFILING
                                    if RANGE_PROFILING
                                        min_inv_r3_term = min([min_inv_r3_term inv_dist_3]);
                                        max_inv_r3_term = max([max_inv_r3_term inv_dist_3]);
                                        min_inv_r6_term = min([min_inv_r6_term inv_dist_6]);
                                        max_inv_r6_term = max([max_inv_r6_term inv_dist_6]);
                                        min_inv_r8_term = min([min_inv_r8_term inv_dist_8]);
                                        max_inv_r8_term = max([max_inv_r8_term inv_dist_8]);
                                        min_inv_r12_term = min([min_inv_r12_term inv_dist_12]);
                                        max_inv_r12_term = max([max_inv_r12_term inv_dist_12]);
                                        min_inv_r14_term = min([min_inv_r14_term inv_dist_14]);
                                        max_inv_r14_term = max([max_inv_r14_term inv_dist_14]);
                                        min_LJ_Force = min([min_LJ_Force LJ_Force_over_R]);
                                        max_LJ_Force = max([max_LJ_Force LJ_Force_over_R]);
                                        min_LJ_Energy = min([min_LJ_Energy LJ_Energy]);
                                        max_LJ_Energy = max([max_LJ_Energy LJ_Energy]);
                                        min_Coulomb_Force = min([min_Coulomb_Force Coulomb_Force_over_R]);
                                        max_Coulomb_Force = max([max_Coulomb_Force Coulomb_Force_over_R]);
                                        min_Coulomb_Energy = min([min_Coulomb_Energy Coulomb_Energy]);
                                        max_Coulomb_Energy = max([max_Coulomb_Energy Coulomb_Energy]);
                                        min_Total_Force = min([min_Total_Force Total_Force_over_R]);
                                        max_Total_Force = max([max_Total_Force Total_Force_over_R]);
                                        min_Total_Energy = min([min_Total_Energy Total_Energy]);
                                        max_Total_Energy = max([max_Total_Energy Total_Energy]);
                                        min_Force_Acc = min([min_Force_Acc Force_Acc_x Force_Acc_y Force_Acc_z]);
                                        max_Force_Acc = max([max_Force_Acc Force_Acc_x Force_Acc_y Force_Acc_z]);
                                        min_Energy_Acc = min([min_Energy_Acc Energy_Acc]);
                                        max_Energy_Acc = max([max_Energy_Acc Energy_Acc]);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                % Assign the force and energy to cell_particle
                cell_particle(home_cell_id,ref_particle_ptr,4) = Force_Acc_x;
                cell_particle(home_cell_id,ref_particle_ptr,5) = Force_Acc_y;
                cell_particle(home_cell_id,ref_particle_ptr,6) = Force_Acc_z;
                cell_particle(home_cell_id,ref_particle_ptr,7) = Energy_Acc;
                cell_particle(home_cell_id,ref_particle_ptr,8) = particles_within_cutoff;
                force_write_back_counter = force_write_back_counter + 1;        % Profiling variable to record how many neighbor particles has been calculated
            end
            
        end
    end
end
fprintf('The number of bonded particle pairs is: %d\n',Bonded_Particle_Pairs/2);
fprintf('Force evaluation is finished!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Profiling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('*** Start Profiling! ***\n')
% Verify the total number of evaluated particles
valid_counter = 0;
max_home_cell_id = 0;
Verification_Bonded_Particle_Pairs = 0;
for home_cell_x = 1:CELL_COUNT_X
    for home_cell_y = 1:CELL_COUNT_Y
        for home_cell_z = 1:CELL_COUNT_Z
            home_cell_id = (home_cell_x-1)*CELL_COUNT_Y*CELL_COUNT_Z + (home_cell_y-1)*CELL_COUNT_Z + home_cell_z;
            if home_cell_id > max_home_cell_id
                max_home_cell_id = home_cell_id;
            end
            particle_count = particle_in_cell_counter(home_cell_x,home_cell_y,home_cell_z);
            for particle_ptr = 1:particle_count
                if cell_particle(home_cell_id,particle_ptr,4) == 0 && cell_particle(home_cell_id,particle_ptr,5) == 0 && cell_particle(home_cell_id,particle_ptr,6) == 0
                    fprintf('home cell %d, particle %d, force value missing!\n', home_cell_id,particle_ptr);
                end
                if cell_particle(home_cell_id,particle_ptr,7) == 0
                    fprintf('home cell %d, particle %d, energy value missing!\n', home_cell_id,particle_ptr);
                end
                if cell_particle(home_cell_id,particle_ptr,4) ~= 0 || cell_particle(home_cell_id,particle_ptr,5) ~= 0 || cell_particle(home_cell_id,particle_ptr,6) ~= 0 || cell_particle(home_cell_id,particle_ptr,7) ~= 0
                    valid_counter = valid_counter + 1;
                end
            end
        end
    end
end
fprintf('Total particle number is %d, particles mapped into cell is %d, force write back counter is %d, parciles with valid force value is %d\n', TOTAL_PARTICLE, TOTAL_PARTICLE-out_range_particle_counter, force_write_back_counter, valid_counter);

if RANGE_PROFILING
    max_range = max([max_dx max_dx_2 max_r2 max_inv_r3_term max_inv_r6_term max_inv_r8_term max_inv_r12_term max_inv_r14_term max_LJ_Force max_LJ_Energy max_Coulomb_Force max_Coulomb_Energy max_Total_Force max_Total_Energy max_Force_Acc max_Energy_Acc]);
    min_range = min([min_dx min_dx_2 min_r2 min_inv_r3_term min_inv_r6_term min_inv_r8_term min_inv_r12_term min_inv_r14_term min_LJ_Force min_LJ_Energy min_Coulomb_Force min_Coulomb_Energy min_Total_Force min_Total_Energy min_Force_Acc min_Energy_Acc]);
    fprintf('The value are ranging from %f to %f.\n', min_range, max_range);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Verification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ENABLE_VERIFICATION
    fprintf('Start Verification, Force model is %s.\n', FORCE_MODEL);
    % Pick a reference paricle and evaluate (!!!!!!!!!!avoid the particles in the corner cells!!!!!!!!)
    home_cell_x = 3;
    home_cell_y = 4;
    home_cell_z = 3;
    particle_id = 2;
    if (home_cell_x == 1 || home_cell_x == CELL_COUNT_X || home_cell_y == 1 || home_cell_y == CELL_COUNT_Y || home_cell_z == 1 || home_cell_z == CELL_COUNT_Z)
        fprintf('Error in verification: Please reselect the reference particle cell (avoid the corner cells)!!!!\n');
        return;
    end
    cell_id = (home_cell_x-1)*CELL_COUNT_Y*CELL_COUNT_Z + (home_cell_y-1)*CELL_COUNT_Z + home_cell_z;
    ref_particle_pos_x = cell_particle(cell_id, particle_id, 1);
    ref_particle_pos_y = cell_particle(cell_id, particle_id, 2);
    ref_particle_pos_z = cell_particle(cell_id, particle_id, 3);
    % Traverse across all the particles in the simulation space
    Force_Acc = single(zeros(1,3));
    particles_within_cutoff = 0;
    Verification_Bonded_Particle_Pairs = uint32(0);
    verification_counter_particle_within_cutoff = zeros(CELL_COUNT_X, CELL_COUNT_Y, CELL_COUNT_Z);
    for i = 1:TOTAL_PARTICLE
        neighbor_particle_pos_x = position_data(i,1);
        neighbor_particle_pos_y = position_data(i,2);
        neighbor_particle_pos_z = position_data(i,3);
        dist_x = single(neighbor_particle_pos_x - ref_particle_pos_x);
        dist_y = single(neighbor_particle_pos_y - ref_particle_pos_y);
        dist_z = single(neighbor_particle_pos_z - ref_particle_pos_z);
        dist_x_2 = dist_x^2;
        dist_y_2 = dist_y^2;
        dist_z_2 = dist_z^2;
        dist_2 = dist_x_2 + dist_y_2 + dist_z_2;

        % Filtering logic and force calculation
        % PROFILING: Count the # of bonded particle pairs
        if dist_2 <= BOND_DISTANCE_2
            Verification_Bonded_Particle_Pairs = Verification_Bonded_Particle_Pairs + 1;
        end
        if dist_2 <= CUTOFF_RADIUS_2 && dist_2 > BOND_DISTANCE_2
            % increment the counter for particles within cutoff radius
            particles_within_cutoff = particles_within_cutoff + 1;
            % check if there are particles within cutoff radius but falling into far cells
            neighbor_cell_x = ceil(neighbor_particle_pos_x / CUTOFF_RADIUS);
            neighbor_cell_y = ceil(neighbor_particle_pos_y / CUTOFF_RADIUS);
            neighbor_cell_z = ceil(neighbor_particle_pos_z / CUTOFF_RADIUS);
            if neighbor_cell_x > home_cell_x + 1 || neighbor_cell_x < home_cell_x - 1 || neighbor_cell_y > home_cell_y + 1 || neighbor_cell_y < home_cell_y - 1 || neighbor_cell_z > home_cell_z + 1 || neighbor_cell_z < home_cell_z - 1
                fprintf('!!!!!Exceptions!! neighbor particle (%f,%f,%f) is not captured by reference particle (%f,%f,%f)....\n',neighbor_particle_pos_x,neighbor_particle_pos_y,neighbor_particle_pos_z,ref_particle_pos_x,ref_particle_pos_y,ref_particle_pos_z);
                return;
            end

%             % recording the # of particles that within the cutoff radius in each cell
%             verification_counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z) = verification_counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z) + 1;
%             if neighbor_cell_x == 3 && neighbor_cell_y == 4 && neighbor_cell_z == 2
%                 temp_verification_cell_particles_within_cutoff(verification_counter_particle_within_cutoff(neighbor_cell_x,neighbor_cell_y,neighbor_cell_z),1:3) = [neighbor_particle_pos_x, neighbor_particle_pos_y, neighbor_particle_pos_z]; 
%             end

            switch FORCE_MODEL
                case 'OpenMM'
                    % Calculate the switch function value (http://docs.openmm.org/6.2.0/userguide/theory.html)
                    dist = sqrt(dist_2);
                    x = (dist - SWITCH_DIST) / (CUTOFF_RADIUS - SWITCH_DIST);
                    if dist >= SWITCH_DIST && dist <= CUTOFF_RADIUS
                        Switch_Val = 1 - 10*x^3 + 15*x^4 - 6*x^5;
                        Switch_Deri = (-30*x^2 + 60*x^3 - 30*x^4) / (CUTOFF_RADIUS - SWITCH_DIST);
                    elseif dist >= 0 && dist < SWITCH_DIST
                        Switch_Val = 1;
                        Switch_Deri = 0;
                    else
                        Switch_Val = 0;
                        Switch_Deri = 0;
                    end
                    inv_dist = 1 / sqrt(dist_2);
                    inv_dist_2 = 1 / dist_2;
                    inv_dist_3 = inv_dist_2 * inv_dist;
                    inv_dist_4 = inv_dist_2 ^ 2;
                    inv_dist_6 = inv_dist_2 * inv_dist_4;
                    inv_dist_8 = inv_dist_4 ^ 2;
                    inv_dist_12 = inv_dist_4 * inv_dist_8;
                    inv_dist_14 = inv_dist_4 * inv_dist_8 * inv_dist_2;
                    sigma_6 = SIGMA^6;
                    sigma_12 = sigma_6^2;
                    chargeProd = ONE_4PI_EPS0 * Q1 * Q2;
                    % Coulomb interaction with cutoff using reaction field approximation
                    krf = INV_CUTOFF_RADIUS_3 * (SOLVENT_DIELECTRIC - 1) / (2 * SOLVENT_DIELECTRIC + 1);
                    crf = INV_CUTOFF_RADIUS * (3 * SOLVENT_DIELECTRIC) / (2 * SOLVENT_DIELECTRIC + 1);
                    % Force (The force calculate here is F/r, for easy calculation of force component on each axis)
                    % LJ force over R
                    LJ_Force_over_R = Switch_Val*4*EPSILON*(12*sigma_12*inv_dist_14 - 6*sigma_6*inv_dist_8);
                    % Apply switch condition for LJ force
                    LJ_Force_over_R = LJ_Force_over_R - Switch_Deri*4*EPSILON*(sigma_12*inv_dist_12 - sigma_6*inv_dist_6)*inv_dist;
                    % Coulomb Force over R
                    Coulomb_Force_over_R = chargeProd * (inv_dist_3 - 2*krf);
                    % Total force over R
                    Total_Force_over_R = LJ_Force_over_R + Coulomb_Force_over_R;
                    % Accumulate the force
                    Force_Acc(1) = Force_Acc(1) + Total_Force_over_R * dist_x;
                    Force_Acc(2) = Force_Acc(2) + Total_Force_over_R * dist_y;
                    Force_Acc(3) = Force_Acc(3) + Total_Force_over_R * dist_z;

                case 'CAAD'
                    % Filtering, and direct force evaluation
                    A = 48 * EPSILON *(SIGMA^12);
                    B = 24 * EPSILON *(SIGMA^6);
                    % Force evaluation
                    inv_dist_2 = 1 / dist_2;
                    inv_dist_4 = inv_dist_2 ^ 2;
                    inv_dist_8 = inv_dist_4 ^ 2;
                    inv_dist_14 = inv_dist_4 * inv_dist_8 * inv_dist_2;
                    % Force calculate (??????????????????? Is the formular here right ??????????????????)
                    Coulomb_Force_over_R = CC * erfc(EWALD_COEF*sqrt(dist_2)) * (1/sqrt(dist_2));
                    LJ_Force_over_R = A * inv_dist_14 - B * inv_dist_8;
                    Total_Force_over_R = Coulomb_Force_over_R + LJ_Force_over_R;

                    % Accumulate the force
                    Force_Acc(1) = Force_Acc(1) + Total_Force_over_R * dist_x;
                    Force_Acc(2) = Force_Acc(2) + Total_Force_over_R * dist_y;
                    Force_Acc(3) = Force_Acc(3) + Total_Force_over_R * dist_z;

                otherwise
                    error('Invalid force model for verification!\n');
            end
        end
    end
    fprintf('The simulated force is (%f,%f,%f) with %d particles in cutoff, total bonded paricle pairs is %d\nThe verification value is (%f,%f,%f) with %d particles in cutoff, bonded particle pairs for this reference particle is %d.\n', cell_particle(cell_id,particle_id,4),cell_particle(cell_id,particle_id,5),cell_particle(cell_id,particle_id,6),cell_particle(cell_id,particle_id,8),Bonded_Particle_Pairs/2, Force_Acc(1),Force_Acc(2),Force_Acc(3),particles_within_cutoff,Verification_Bonded_Particle_Pairs);
end