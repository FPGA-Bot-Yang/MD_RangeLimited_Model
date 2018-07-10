%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate LJ Interpolation table for MiniMD code for STX
% Targeting the force_compute function in MiniMD OpenCL code
% There's no index of A, B used, the force model is F_over_r = 48*R^-14 + 24*R^-8
% Accumuation is based on the component in each direction: F_acc_x = F_over_r*dx, F_acc_y = F_over_r*dy, F_acc_z = F_over_r*dz
%
%% Code reference
% __kernel void force_compute(_constant MMD_floatK3* _restrict  x, __global MMD_floatK3* _restrict f, __global int* _restrict numneigh,
%                                         __global int* _restrict neighbors, int maxneighs, int nlocal, MMD_float cutforcesq)
% {
%  for(int i=0;i<nlocal;i++)
%  {
% 
%    __global int* neighs = neighbors + i;
%    MMD_floatK3 ftmp;
%    MMD_floatK3 xi = x[i];
%    printf("I=%d;(%d,%d,%d)\n",i,xi.x,xi.y,xi.z);
%    MMD_floatK3 fi = {0.0f,0.0f,0.0f};
% 
%    for (int k = 0; k < numneigh[i]; k++) {
%      int j = neighs[k*nlocal];
%      MMD_floatK3 delx = xi - x[j];
%      MMD_float rsq = delx.x*delx.x + delx.y*delx.y + delx.z*delx.z;
%      if (rsq < cutforcesq) {
%       MMD_float sr2 = 1.0f/rsq;
%       MMD_float sr6 = sr2*sr2*sr2;
%       MMD_float force = 48.0f*sr6*(sr6-0.5f)*sr2;
%       fi += force * delx;
%      }
%    }
%    f[i] = fi;
%  }
% }
%
% By: Chen Yang
% 07/09/2018
% Boston University, CAAD Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example:
% Compute coefficients for target function f=x^-7 on range from 2^0 to 2^8 for a 2nd or 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 or 256 intervals, and output results to a.txt.
% 2^0 is becuase r2 > exclude_r2_min
% 2^8 is becuase r2 < cutoff2
%
% Command line:
% interpolation_order = 1;                % interpolation order, no larger than 3
% segment_num = 12;                       % # of segment
% bin_num = 256;                          % # of bins per segment
% precision = 8;                          % # of datepoints when generating the polynomial index 
% min_range = 1;                          % minimal range for the evaluation
% max_range = min_range * 2^segment_num;  % maximum range for the evaluation (currently this is the cutoff radius)
% MiniMD_LJ_poly_interpolation_function(3,12,256,8,1,2^12)


% interpolation_order: interpolation order
% bin_num: # of bins per segment
% precision: # of datapoints for each interpolation
% min, max: range of distance


function  MiniMD_LJ_poly_interpolation_function(interpolation_order,segment_num,bin_num,precision,min,max)
% interpolation_order is the order of interpolation. i.e, interpolation_order=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if nargin ~= 6
        error('LJ_poly_interpolation_function(interpolation_order,segment_num,bin_num,precision,min,max)');
    end
    
	if min < 0
		error('min must be greater than 0.');
    end
    
	if interpolation_order > 3 || interpolation_order <= 0
        error('The supported interpolation order is 1, 2, 3 ....\n');
	end

    %% Create output files
    % Force term
%     fileID_output_force_14 = fopen('force_value_r14.txt', 'wt');
%     fileID_output_force_8 = fopen('force_value_r8.txt', 'wt');
    
	fileID_0  = fopen('table_c0_r14.txt', 'wt');
    fileID_1  = fopen('table_c1_r14.txt', 'wt');
    fileID_mif_0 = fopen('table_c0_r14.mif', 'wt');
    fileID_mif_1 = fopen('table_c1_r14.mif', 'wt');
    if(interpolation_order > 1)
        fileID_2  = fopen('table_c2_r14.txt', 'wt');
        fileID_mif_2 = fopen('table_c2_r14.mif', 'wt');
    end
    if(interpolation_order > 2)
        fileID_3  = fopen('table_c3_r14.txt', 'wt');
        fileID_mif_3 = fopen('table_c3_r14.mif', 'wt');
    end
    
    fileID_4  = fopen('table_c0_r8.txt', 'wt');
    fileID_5  = fopen('table_c1_r8.txt', 'wt');
    fileID_mif_4 = fopen('table_c0_r8.mif', 'wt');
    fileID_mif_5 = fopen('table_c1_r8.mif', 'wt');
    if(interpolation_order > 1)
        fileID_6  = fopen('table_c2_r8.txt', 'wt');
        fileID_mif_6 = fopen('table_c2_r8.mif', 'wt');
    end
    if(interpolation_order > 2)
        fileID_7  = fopen('table_c3_r8.txt', 'wt');
        fileID_mif_7 = fopen('table_c3_r8.mif', 'wt');
    end
    
    % Energy term
    fileID_8  = fopen('table_c0_r12.txt', 'wt');
    fileID_9  = fopen('table_c1_r12.txt', 'wt');
    fileID_mif_8 = fopen('table_c0_r12.mif', 'wt');
    fileID_mif_9 = fopen('table_c1_r12.mif', 'wt');
    if(interpolation_order > 1)
        fileID_10 = fopen('table_c2_r12.txt', 'wt');
        fileID_mif_10 = fopen('table_c2_r12.mif', 'wt');
    end
    if(interpolation_order > 2)
        fileID_11 = fopen('table_c3_r12.txt', 'wt');
        fileID_mif_11 = fopen('table_c3_r12.mif', 'wt');
    end
    
    fileID_12 = fopen('table_c0_r6.txt', 'wt');
    fileID_13 = fopen('table_c1_r6.txt', 'wt');
    fileID_mif_12 = fopen('table_c0_r6.mif', 'wt');
    fileID_mif_13 = fopen('table_c1_r6.mif', 'wt');
    if(interpolation_order > 1)
        fileID_14 = fopen('table_c2_r6.txt', 'wt');
        fileID_mif_14 = fopen('table_c2_r6.mif', 'wt');
    end
    if(interpolation_order > 2)
        fileID_15 = fopen('table_c3_r6.txt', 'wt');
        fileID_mif_15 = fopen('table_c3_r6.mif', 'wt');
    end

    %% Write the mif file header
    fprintf(fileID_mif_0,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_0,'WIDTH = 32;\n');
    fprintf(fileID_mif_0,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_0,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_0,'CONTENT\n');
    fprintf(fileID_mif_0,'BEGIN\n');
    
    fprintf(fileID_mif_1,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_1,'WIDTH = 32;\n');
    fprintf(fileID_mif_1,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_1,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_1,'CONTENT\n');
    fprintf(fileID_mif_1,'BEGIN\n');
    
    fprintf(fileID_mif_4,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_4,'WIDTH = 32;\n');
    fprintf(fileID_mif_4,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_4,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_4,'CONTENT\n');
    fprintf(fileID_mif_4,'BEGIN\n');
    
    fprintf(fileID_mif_5,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_5,'WIDTH = 32;\n');
    fprintf(fileID_mif_5,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_5,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_5,'CONTENT\n');
    fprintf(fileID_mif_5,'BEGIN\n');
    
    fprintf(fileID_mif_8,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_8,'WIDTH = 32;\n');
    fprintf(fileID_mif_8,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_8,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_8,'CONTENT\n');
    fprintf(fileID_mif_8,'BEGIN\n');
    
    fprintf(fileID_mif_9,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_9,'WIDTH = 32;\n');
    fprintf(fileID_mif_9,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_9,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_9,'CONTENT\n');
    fprintf(fileID_mif_9,'BEGIN\n');
    
    fprintf(fileID_mif_12,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_12,'WIDTH = 32;\n');
    fprintf(fileID_mif_12,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_12,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_12,'CONTENT\n');
    fprintf(fileID_mif_12,'BEGIN\n');
    
    fprintf(fileID_mif_13,'DEPTH = %d;\n',segment_num*bin_num);
    fprintf(fileID_mif_13,'WIDTH = 32;\n');
    fprintf(fileID_mif_13,'ADDRESS_RADIX = DEC;\n');
    fprintf(fileID_mif_13,'DATA_RADIX = HEX;\n');
    fprintf(fileID_mif_13,'CONTENT\n');
    fprintf(fileID_mif_13,'BEGIN\n');
    
    if interpolation_order > 1
        fprintf(fileID_mif_2,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_2,'WIDTH = 32;\n');
        fprintf(fileID_mif_2,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_2,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_2,'CONTENT\n');
        fprintf(fileID_mif_2,'BEGIN\n');

        fprintf(fileID_mif_6,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_6,'WIDTH = 32;\n');
        fprintf(fileID_mif_6,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_6,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_6,'CONTENT\n');
        fprintf(fileID_mif_6,'BEGIN\n');
        
        fprintf(fileID_mif_10,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_10,'WIDTH = 32;\n');
        fprintf(fileID_mif_10,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_10,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_10,'CONTENT\n');
        fprintf(fileID_mif_10,'BEGIN\n');

        fprintf(fileID_mif_14,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_14,'WIDTH = 32;\n');
        fprintf(fileID_mif_14,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_14,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_14,'CONTENT\n');
        fprintf(fileID_mif_14,'BEGIN\n');
    end
    
    if interpolation_order > 2
        fprintf(fileID_mif_3,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_3,'WIDTH = 32;\n');
        fprintf(fileID_mif_3,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_3,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_3,'CONTENT\n');
        fprintf(fileID_mif_3,'BEGIN\n');

        fprintf(fileID_mif_7,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_7,'WIDTH = 32;\n');
        fprintf(fileID_mif_7,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_7,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_7,'CONTENT\n');
        fprintf(fileID_mif_7,'BEGIN\n');
        
        fprintf(fileID_mif_11,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_11,'WIDTH = 32;\n');
        fprintf(fileID_mif_11,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_11,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_11,'CONTENT\n');
        fprintf(fileID_mif_11,'BEGIN\n');

        fprintf(fileID_mif_15,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_15,'WIDTH = 32;\n');
        fprintf(fileID_mif_15,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_15,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_15,'CONTENT\n');
        fprintf(fileID_mif_15,'BEGIN\n');
    end
    
    
    %% Start evaluation
    count = 0;
    range_min = min;
    range_max = 2*min;
    while(range_min < max)          % EACH SEGMENT
    
        step = single((range_max-range_min)/bin_num);
        ca = single(range_min);
    
        for j=1:bin_num                   % EACH BIN
            x = zeros(precision,1);

            % for L-J Potential/Force
            inv_r12 = single(zeros(precision,1));
            inv_r6  = single(zeros(precision,1));
            inv_r14 = single(zeros(precision,1));
            inv_r8  = single(zeros(precision,1));
            final_inv_r6 = single(zeros(precision,1));
            final_inv_r12 = single(zeros(precision,1));
            final_inv_r8 = single(zeros(precision,1));
            final_inv_r14 = single(zeros(precision,1));
            
            delta = step/precision;
    
            for i=1:precision
                x(i) = ca + (i-1) * delta;
                r2 = x(i);

                inv_r2 = 1/r2;
                inv_r4 = inv_r2 * inv_r2;
                
                inv_r6(i)  = inv_r2 * inv_r4;
                inv_r12(i) = inv_r6(i) * inv_r6(i);
                inv_r14(i) = inv_r12(i) * inv_r2;
                inv_r8(i) = inv_r6(i) * inv_r2;
                
                % final term with coefficients
                final_inv_r6(i) = 4 * inv_r6(i);
                final_inv_r12(i) = 4 * inv_r12(i);
                final_inv_r14(i) = 48 * inv_r14(i);
                final_inv_r8(i)  = 24 * inv_r8(i);
                
%                 fprintf(fileID_output_force_14,'%15.25f\n',final_inv_r14(i));
%                 fprintf(fileID_output_force_8,'%15.25f\n',final_inv_r8(i));
            end

            r14_func = polyfit(x,final_inv_r14,interpolation_order);
            r8_func  = polyfit(x,final_inv_r8,interpolation_order);
            r12_func = polyfit(x,final_inv_r12,interpolation_order);
            r6_func  = polyfit(x,final_inv_r6,interpolation_order);
            
            switch(interpolation_order)
                case 1
                    % write to file for verification
                    fprintf(fileID_0,'%15.25f\n',r14_func(2));
                    fprintf(fileID_1,'%15.25f\n',r14_func(1));

                    fprintf(fileID_4,'%15.25f\n',r8_func(2));
                    fprintf(fileID_5,'%15.25f\n',r8_func(1));

                    fprintf(fileID_8, '%15.25f\n',r12_func(2));
                    fprintf(fileID_9, '%15.25f\n',r12_func(1));

                    fprintf(fileID_12,'%15.25f\n',r6_func(2));
                    fprintf(fileID_13,'%15.25f\n',r6_func(1));
                    
                    % write to file for mif
                    fprintf(fileID_mif_0,'%d : %tX;\n',count, r14_func(2));
                    fprintf(fileID_mif_1,'%d : %tX;\n',count, r14_func(1));

                    fprintf(fileID_mif_4,'%d : %tX;\n',count, r8_func(2));
                    fprintf(fileID_mif_5,'%d : %tX;\n',count, r8_func(1));

                    fprintf(fileID_mif_8,'%d : %tX;\n',count, r12_func(2));
                    fprintf(fileID_mif_9,'%d : %tX;\n',count, r12_func(1));

                    fprintf(fileID_mif_12,'%d : %tX;\n',count, r6_func(2));
                    fprintf(fileID_mif_13,'%d : %tX;\n',count, r6_func(1));
                    
                case 2
                    % write to file for verification
                    fprintf(fileID_0,'%15.25f\n',r14_func(3));
                    fprintf(fileID_1,'%15.25f\n',r14_func(2));
                    fprintf(fileID_2,'%15.25f\n',r14_func(1)); 

                    fprintf(fileID_4,'%15.25f\n',r8_func(3));
                    fprintf(fileID_5,'%15.25f\n',r8_func(2));
                    fprintf(fileID_6,'%15.25f\n',r8_func(1));

                    fprintf(fileID_8, '%15.25f\n',r12_func(3));
                    fprintf(fileID_9, '%15.25f\n',r12_func(2));
                    fprintf(fileID_10,'%15.25f\n',r12_func(1)); 

                    fprintf(fileID_12,'%15.25f\n',r6_func(3));
                    fprintf(fileID_13,'%15.25f\n',r6_func(2));
                    fprintf(fileID_14,'%15.25f\n',r6_func(1));
                    
                    % write to file for mif
                    fprintf(fileID_mif_0,'%d : %tX;\n',count, r14_func(3));
                    fprintf(fileID_mif_1,'%d : %tX;\n',count, r14_func(2));
                    fprintf(fileID_mif_2,'%d : %tX;\n',count, r14_func(1));

                    fprintf(fileID_mif_4,'%d : %tX;\n',count, r8_func(3));
                    fprintf(fileID_mif_5,'%d : %tX;\n',count, r8_func(2));
                    fprintf(fileID_mif_6,'%d : %tX;\n',count, r8_func(1));

                    fprintf(fileID_mif_8,'%d : %tX;\n',count, r12_func(3));
                    fprintf(fileID_mif_9,'%d : %tX;\n',count, r12_func(2));
                    fprintf(fileID_mif_10,'%d : %tX;\n',count, r12_func(1));

                    fprintf(fileID_mif_12,'%d : %tX;\n',count, r6_func(3));
                    fprintf(fileID_mif_13,'%d : %tX;\n',count, r6_func(2));
                    fprintf(fileID_mif_14,'%d : %tX;\n',count, r6_func(1));
                    
                case 3
                    % write to file for verification
                    fprintf(fileID_0,'%15.25f\n',r14_func(4));
                    fprintf(fileID_1,'%15.25f\n',r14_func(3));
                    fprintf(fileID_2,'%15.25f\n',r14_func(2)); 
                    fprintf(fileID_3,'%15.25f\n',r14_func(1)); 

                    fprintf(fileID_4,'%15.25f\n',r8_func(4));
                    fprintf(fileID_5,'%15.25f\n',r8_func(3));
                    fprintf(fileID_6,'%15.25f\n',r8_func(2));
                    fprintf(fileID_7,'%15.25f\n',r8_func(1));

                    fprintf(fileID_8, '%15.25f\n',r12_func(4));
                    fprintf(fileID_9, '%15.25f\n',r12_func(3));
                    fprintf(fileID_10,'%15.25f\n',r12_func(2)); 
                    fprintf(fileID_11,'%15.25f\n',r12_func(1));

                    fprintf(fileID_12,'%15.25f\n',r6_func(4));
                    fprintf(fileID_13,'%15.25f\n',r6_func(3));
                    fprintf(fileID_14,'%15.25f\n',r6_func(2));
                    fprintf(fileID_15,'%15.25f\n',r6_func(1));
                    
                    % write to file for mif
                    fprintf(fileID_mif_0,'%d : %tX;\n',count, r14_func(4));
                    fprintf(fileID_mif_1,'%d : %tX;\n',count, r14_func(3));
                    fprintf(fileID_mif_2,'%d : %tX;\n',count, r14_func(2));
                    fprintf(fileID_mif_3,'%d : %tX;\n',count, r14_func(1));

                    fprintf(fileID_mif_4,'%d : %tX;\n',count, r8_func(4));
                    fprintf(fileID_mif_5,'%d : %tX;\n',count, r8_func(3));
                    fprintf(fileID_mif_6,'%d : %tX;\n',count, r8_func(2));
                    fprintf(fileID_mif_7,'%d : %tX;\n',count, r8_func(1));

                    fprintf(fileID_mif_8,'%d : %tX;\n',count, r12_func(4));
                    fprintf(fileID_mif_9,'%d : %tX;\n',count, r12_func(3));
                    fprintf(fileID_mif_10,'%d : %tX;\n',count, r12_func(2));
                    fprintf(fileID_mif_11,'%d : %tX;\n',count, r12_func(1));

                    fprintf(fileID_mif_12,'%d : %tX;\n',count, r6_func(4));
                    fprintf(fileID_mif_13,'%d : %tX;\n',count, r6_func(3));
                    fprintf(fileID_mif_14,'%d : %tX;\n',count, r6_func(2));
                    fprintf(fileID_mif_15,'%d : %tX;\n',count, r6_func(1));
            end
            
            %fprintf(fileID_0,'%d\t:\t%tx;\n', count, r14_func(4));
            %fprintf(fileID_1,'%d\t:\t%tx;\n', count, r14_func(3));
            %fprintf(fileID_2,'%d\t:\t%tx;\n', count, r14_func(2));
            %fprintf(fileID_3,'%d\t:\t%tx;\n', count, r14_func(1));
            
            %fprintf(fileID_4,'%d\t:\t%tx;\n', count, r8_func(4));
            %fprintf(fileID_5,'%d\t:\t%tx;\n', count, r8_func(3));
            %fprintf(fileID_6,'%d\t:\t%tx;\n', count, r8_func(2));
            %fprintf(fileID_7,'%d\t:\t%tx;\n', count, r8_func(1));
            
            %fprintf(fileID_8,'%d\t:\t%tx;\n',  count, r12_func(4));
            %fprintf(fileID_9,'%d\t:\t%tx;\n',  count, r12_func(3));
            %fprintf(fileID_10,'%d\t:\t%tx;\n', count, r12_func(2));
            %fprintf(fileID_11,'%d\t:\t%tx;\n', count, r12_func(1));
            
            %fprintf(fileID_12,'%d\t:\t%tx;\n', count, r6_func(4));
            %fprintf(fileID_13,'%d\t:\t%tx;\n', count, r6_func(3));
            %fprintf(fileID_14,'%d\t:\t%tx;\n', count, r6_func(2));
            %fprintf(fileID_15,'%d\t:\t%tx;\n', count, r6_func(1));
           
            count = count + 1;
            ca = ca + step;
        end
        
        range_min = range_min * 2;
        range_max = range_max * 2;
    end
    
    %% Write the end of mif file
    fprintf(fileID_mif_0,'END;\n');
    fprintf(fileID_mif_1,'END;\n');
    fprintf(fileID_mif_4,'END;\n');
    fprintf(fileID_mif_5,'END;\n');
    if interpolation_order > 1
        fprintf(fileID_mif_2,'END;\n');
        fprintf(fileID_mif_6,'END;\n');
    end
    if interpolation_order > 2
        fprintf(fileID_mif_3,'END;\n');
        fprintf(fileID_mif_7,'END;\n');
    end
    
    
    %% Close the files
%     fclose(fileID_output_force_14);
%     fclose(fileID_output_force_8);
    
    fclose(fileID_0);
    fclose(fileID_1);
    fclose(fileID_4);
    fclose(fileID_5);
    fclose(fileID_8);
    fclose(fileID_9);
    fclose(fileID_12);
    fclose(fileID_13);
    
    fclose(fileID_mif_0);
    fclose(fileID_mif_1);
    fclose(fileID_mif_4);
    fclose(fileID_mif_5);
    fclose(fileID_mif_8);
    fclose(fileID_mif_9);
    fclose(fileID_mif_12);
    fclose(fileID_mif_13);

    if interpolation_order > 1
        fclose(fileID_2);
        fclose(fileID_6);
        fclose(fileID_10);
        fclose(fileID_14);
        
        fclose(fileID_mif_2);
        fclose(fileID_mif_6);
        fclose(fileID_mif_10);
        fclose(fileID_mif_14);
    end
    
    if interpolation_order > 2
        fclose(fileID_3);
        fclose(fileID_7);
        fclose(fileID_11);
        fclose(fileID_15);
        
        fclose(fileID_mif_3);
        fclose(fileID_mif_7);
        fclose(fileID_mif_11);
        fclose(fileID_mif_15);
    end
    