% Example:
% Compute coefficients for target function f=x^-7 on range from 2^0 to 2^8 for a 2nd or 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 or 256 intervals, and output results to a.txt.
% 2^0 is becuase r2 > exclude_r2_min
% 2^8 is becuase r2 < cutoff2
%
% Command line:
% syms x
% f=x^-7
% poly_interpolation_coef_real(3,16,2^-4,2^-3,'a.txt')


% m: interpolation order
% n: # of bins per segment
% precision: # of datapoints for each interpolation
% min, max: range of distance
% cutoff: cut-off radius
% switchon: switch on distance
function  LJ_poly_interpolation_function(m,n,precision,min,max,cutoff,switchon)
% m is the order of interpolation. i.e, m=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if(nargin<5)
		error('poly_interpolation_function(m,n,precision,min,max)');
	end
	if(min<0)
		error('min must be greater than 0.');
	end

	fileID_0  = fopen('file_c0_vdw14_f_order3_bin64.txt', 'wt');
    fileID_1  = fopen('file_c1_vdw14_f_order3_bin64.txt', 'wt');
    if(m > 1)
        fileID_2  = fopen('file_c2_vdw14_f_order3_bin64.txt', 'wt');
    end
    if(m > 2)
        fileID_3  = fopen('file_c3_vdw14_f_order3_bin64.txt', 'wt');
    end
    
    fileID_4  = fopen('file_c0_vdw8_f_order3_bin64.txt', 'wt');
    fileID_5  = fopen('file_c1_vdw8_f_order3_bin64.txt', 'wt');
    if(m > 1)
        fileID_6  = fopen('file_c2_vdw8_f_order3_bin64.txt', 'wt');
    end
    if(m > 2)
        fileID_7  = fopen('file_c3_vdw8_f_order3_bin64.txt', 'wt');
    end
    
    fileID_8  = fopen('file_c0_vdw12_e_order3_bin64.txt', 'wt');
    fileID_9  = fopen('file_c1_vdw12_e_order3_bin64.txt', 'wt');
    if(m > 1)
        fileID_10 = fopen('file_c2_vdw12_e_order3_bin64.txt', 'wt');
    end
    if(m > 2)
        fileID_11 = fopen('file_c3_vdw12_e_order3_bin64.txt', 'wt');
    end
    
    fileID_12 = fopen('file_c0_vdw6_e_order3_bin64.txt', 'wt');
    fileID_13 = fopen('file_c1_vdw6_e_order3_bin64.txt', 'wt');
    if(m > 1)
        fileID_14 = fopen('file_c2_vdw6_e_order3_bin64.txt', 'wt');
    end
    if(m > 2)
        fileID_15 = fopen('file_c3_vdw6_e_order3_bin64.txt', 'wt');
    end
    
    count = 0;
    
    cutoff2 = single(cutoff * cutoff);
    switchon2 = single(switchon * switchon);
    inv_denom = single((cutoff2 - switchon2)^3);
    denom = 1/inv_denom;
    
    range_min = min;
    range_max = 2*min;
    
    while(range_min < max)          % EACH SEGMENT
    
        step = single((range_max-range_min)/n);
        ca = single(range_min);
        cb = range_min + step;
    
        for j=1:n                   % EACH BIN
            x = zeros(precision,1);

            % for L-J Potential/Force
            inv_r12 = single(zeros(precision,1));
            inv_r6  = single(zeros(precision,1));
            
            inv_r14 = single(zeros(precision,1));
            inv_r8  = single(zeros(precision,1));
            
            delta = step/precision;
    
            for i=1:precision
                x(i) = ca + (i-1) * delta;
                r2 = x(i);
                
                if(r2 <= switchon2)
                    s = 1;
                    ds = 0;
                end
        
                if(r2 > switchon2 && r2 <= cutoff2)
                    s = (cutoff2 - r2) * (cutoff2 - r2) * (cutoff2 + 2*r2 - 3 * switchon2) * denom;
                    ds = 12 * (cutoff2 - r2) * (switchon2 - r2) * denom;
                end
        
                if(r2 > cutoff2)
                    s = 0;
                    ds = 0;
                end
        
                inv_r2 = 1/r2;
                inv_r4 = inv_r2 * inv_r2;
                
                inv_r6(i)  = inv_r2 * inv_r4 * s;
                inv_r12(i) = inv_r6(i) * inv_r6(i) * s;
                
                inv_r14(i) = inv_r12(i) * (ds - 12*s*inv_r2);
                inv_r8(i)  = inv_r6(i)  * (ds -  6*s*inv_r2);
            end

            r14_func = polyfit(x,inv_r14,m);
            r8_func  = polyfit(x,inv_r8,m);
            r12_func = polyfit(x,inv_r12,m);
            r6_func  = polyfit(x,inv_r6,m);
            
            ncoef=length(r14_func);
            
            switch(m)
                case 1
                    fprintf(fileID_0,'%15.25f\n',r14_func(2));
                    fprintf(fileID_1,'%15.25f\n',r14_func(1));

                    fprintf(fileID_4,'%15.25f\n',r8_func(2));
                    fprintf(fileID_5,'%15.25f\n',r8_func(1));

                    fprintf(fileID_8, '%15.25f\n',r12_func(2));
                    fprintf(fileID_9, '%15.25f\n',r12_func(1));

                    fprintf(fileID_12,'%15.25f\n',r6_func(2));
                    fprintf(fileID_13,'%15.25f\n',r6_func(1));
                    
                case 2
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
                    
                case 3
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
    
    fclose(fileID_0);
    fclose(fileID_1);
    fclose(fileID_4);
    fclose(fileID_5);
    fclose(fileID_8);
    fclose(fileID_9);
    fclose(fileID_12);
    fclose(fileID_13);

    if m > 1
        fclose(fileID_2);
        fclose(fileID_6);
        fclose(fileID_10);
        fclose(fileID_14);
    end
    
    if m > 2
        fclose(fileID_3);
        fclose(fileID_7);
        fclose(fileID_11);
        fclose(fileID_15);
    end
    