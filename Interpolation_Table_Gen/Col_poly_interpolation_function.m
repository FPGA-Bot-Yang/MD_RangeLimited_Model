% Example:
% Compute coefficients for target functio` f=x^-7 on range from 2^0 to 2^8 for a 2nd or 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 or 256 intervals, and output results to a.txt.
% 2^0 is becuase r2 > exclude_r2_min
% 2^8 is becuase r2 < cutoff2
%
% Command line:
% syms x
% f=x^-7
% poly_interpolation_coef_real(3,16,2^-4,2^-3,'a.txt')

function  Col_poly_interpolation_function(interpolation_order,segment_num,bin_num,precision,min,max)
% interpolation_order is the order of interpolation. i.e, interpolation_order=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
% bin_num is the number of small sections in each segment we have
	if nargin < 6
        error('Col_poly_interpolation_function(interpolation_order,segment_num,bin_num,precision,min,max)');
	end
	if min < 0
        error('min must be gt 0.');
	end
	if interpolation_order > 3 || interpolation_order <= 0
        error('The supported interpolation order is 1, 2, 3 ....\n');
	end

    %% Create output files
	fileID_0 = fopen('file_c0_elec_f.txt', 'wt');
    fileID_1 = fopen('file_c1_elec_f.txt', 'wt');
    fileID_mif_0 = fopen('file_c0_elec_f.mif', 'wt');
    fileID_mif_1 = fopen('file_c1_elec_f.mif', 'wt');
    if interpolation_order > 1
        fileID_2 = fopen('file_c2_elec_f.txt', 'wt');
        fileID_mif_2 = fopen('file_c2_elec_f.mif', 'wt');
    end
    if interpolation_order > 2
        fileID_3 = fopen('file_c3_elec_f.txt', 'wt');
        fileID_mif_3 = fopen('file_c3_elec_f.mif', 'wt');
    end
    
    fileID_4 = fopen('file_c0_elec_e.txt', 'wt');
    fileID_5 = fopen('file_c1_elec_e.txt', 'wt');
    if interpolation_order > 1
        fileID_6 = fopen('file_c2_elec_e.txt', 'wt');
    end
    if interpolation_order > 2
        fileID_7 = fopen('file_c3_elec_e.txt', 'wt');
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
    
    if interpolation_order > 1
        fprintf(fileID_mif_2,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_2,'WIDTH = 32;\n');
        fprintf(fileID_mif_2,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_2,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_2,'CONTENT\n');
        fprintf(fileID_mif_2,'BEGIN\n');
    end
    
    if interpolation_order > 2
        fprintf(fileID_mif_3,'DEPTH = %d;\n',segment_num*bin_num);
        fprintf(fileID_mif_3,'WIDTH = 32;\n');
        fprintf(fileID_mif_3,'ADDRESS_RADIX = DEC;\n');
        fprintf(fileID_mif_3,'DATA_RADIX = HEX;\n');
        fprintf(fileID_mif_3,'CONTENT\n');
        fprintf(fileID_mif_3,'BEGIN\n');
    end
    
    
    %% Start evaluation
    count = 0;
    
    ewald_coef = 0.257952;
    ewald_coef2 = ewald_coef * ewald_coef;
    cc = 332.0636;
    grad_coef = 0.291067;
    
    range_min = min;
    range_max = 2*min;
    
    while(range_min < max)
    
        step = (range_max-range_min)/bin_num;
        ca = range_min;
        %cb = range_min + step;
    
        for j=1:bin_num
            x = zeros(precision,1);
            inv_r3 = zeros(precision,1);
            inv_r  = zeros(precision,1); 
            
            delta = step/precision;
    
            for i=1:precision
                x(i) = ca + (i-1) * delta;
                k = x(i);
                  
                %inv_r3(i) = -1 * cc * (1/k) *( erfc(ewald_coef*sqrt(k)) * (1/sqrt(k)) + grad_coef * exp(-1*k*ewald_coef2) );
                inv_r(i) = cc * erfc(ewald_coef*sqrt(k)) * (1/sqrt(k));
                inv_r3(i) = -1 * (1/k) * ( cc * grad_coef * exp(-1*k*ewald_coef2) + inv_r(i) );     
            end
    
            r3_func = polyfit(x,inv_r3,interpolation_order);
            r_func = polyfit(x,inv_r,interpolation_order);
            
            ncoef=length(r3_func);
            
            switch(interpolation_order)
                case 1
                    fprintf(fileID_0,'%15.25f\n', r3_func(2));
                    fprintf(fileID_1,'%15.25f\n', r3_func(1));

                    fprintf(fileID_4,'%15.25f\n', r_func(2));
                    fprintf(fileID_5,'%15.25f\n', r_func(1));

                    fprintf(fileID_mif_0,'%d : %tX;\n', count, r3_func(2));
                    fprintf(fileID_mif_1,'%d : %tX;\n', count, r3_func(1));

                    
                case 2
                    fprintf(fileID_0,'%15.25f\n', r3_func(3));
                    fprintf(fileID_1,'%15.25f\n', r3_func(2));
                    fprintf(fileID_2,'%15.25f\n', r3_func(1));

                    fprintf(fileID_4,'%15.25f\n', r_func(3));
                    fprintf(fileID_5,'%15.25f\n', r_func(2));
                    fprintf(fileID_6,'%15.25f\n', r_func(1));

                    fprintf(fileID_mif_0,'%d : %tX;\n', count, r3_func(3));
                    fprintf(fileID_mif_1,'%d : %tX;\n', count, r3_func(2));
                    fprintf(fileID_mif_2,'%d : %tX;\n', count, r3_func(1));
                    
                case 3
                    fprintf(fileID_0,'%15.25f\n', r3_func(4));
                    fprintf(fileID_1,'%15.25f\n', r3_func(3));
                    fprintf(fileID_2,'%15.25f\n', r3_func(2));
                    fprintf(fileID_3,'%15.25f\n', r3_func(1));

                    fprintf(fileID_4,'%15.25f\n', r_func(4));
                    fprintf(fileID_5,'%15.25f\n', r_func(3));
                    fprintf(fileID_6,'%15.25f\n', r_func(2));
                    fprintf(fileID_7,'%15.25f\n', r_func(1));

                    fprintf(fileID_mif_0,'%d : %tX;\n', count, r3_func(4));
                    fprintf(fileID_mif_1,'%d : %tX;\n', count, r3_func(3));
                    fprintf(fileID_mif_2,'%d : %tX;\n', count, r3_func(2));
                    fprintf(fileID_mif_3,'%d : %tX;\n', count, r3_func(1));
            end
            
            
            %fprintf(fileID_0,'%d\t:\t%tx;\n', count, r3_func(3));
            %fprintf(fileID_1,'%d\t:\t%tx;\n', count, r3_func(2));
            %fprintf(fileID_2,'%d\t:\t%tx;\n', count, r3_func(1));
            
            %fprintf(fileID_3,'%d\t:\t%tx;\n', count, r_func(3));
            %fprintf(fileID_4,'%d\t:\t%tx;\n', count, r_func(2));
            %fprintf(fileID_5,'%d\t:\t%tx;\n', count, r_func(1));
             
            count = count + 1;
            ca = ca + step;
        end
        range_min = range_min * 2;
        range_max = range_max * 2;
    end
    
    %% Write the end of mif file
    fprintf(fileID_mif_0,'END;\n');
    fprintf(fileID_mif_1,'END;\n');
    if interpolation_order > 1
        fprintf(fileID_mif_2,'END;\n');
    end
    if interpolation_order > 2
        fprintf(fileID_mif_3,'END;\n');
    end
    
    
    %% Close the files
    fclose(fileID_0);
    fclose(fileID_1);
    fclose(fileID_4);
    fclose(fileID_5);
    fclose(fileID_mif_0);
    fclose(fileID_mif_1);
    
    if interpolation_order > 1
        fclose(fileID_2);
        fclose(fileID_6);
        fclose(fileID_mif_2);
    end
    if interpolation_order > 1
        fclose(fileID_3);
        fclose(fileID_7);
        fclose(fileID_mif_3);
    end
    
    
    
   
    