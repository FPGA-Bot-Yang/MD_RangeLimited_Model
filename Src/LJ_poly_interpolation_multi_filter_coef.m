% Example:
% Compute coefficients for target function f=x^-7 on range from 2^0 to 2^8 for a 2nd or 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 or 256 intervals, and output results to a.txt.
%
% Command line:
% syms x
% f=x^-7
% LJ_poly_interpolation_multi_filter_coef(3,16,2^-4,2^-3,'a.txt')

function  LJ_poly_interpolation_multi_filter_coef(m,n,min,max)
% m is the order of interpolation. i.e, m=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if(nargin<4)
		error('LJ_poly_interpolation_multi_filter_coef(m,n,min,max)');
	end
	if(min<0)
		error('min must be gt 0.');
	end

	fileID1 = fopen('s_coef.txt', 'wt');
    fileID2 = fopen('ds_coef.txt', 'wt');
    fileID3 = fopen('r6_order0_coef0.txt', 'wt');
    fileID4 = fopen('r12_order0_coef0.txt', 'wt');
    fileID5 = fopen('r8_order0_coef0.txt', 'wt');
    fileID6 = fopen('r14_order0_coef0.txt', 'wt');
    fileID7 = fopen('r_order0_coef0.txt', 'wt');
    fileID8 = fopen('r3_order0_coef0.txt', 'wt');
    
    % VDW coefficients
    cutoff = 12;
    switchon = 10;
    
    cutoff2 = cutoff * cutoff;
    switchon2 = switchon * switchon;
    inv_denom = (cutoff2 - switchon2)^3;
    denom = 1/inv_denom;
    count = 0;
    
    % PME coefficients
    ewald_coef = 0.257952;
    ewald_coef2 = ewald_coef * ewald_coef;
    cc = 332.0636;
    grad_coef = 0.291067;
    
    range_min = min;
    range_max = 2*min;
    
    while(range_min < max)
        
        x = zeros(n,1);
    
        delta = (range_max-range_min)/n;
       
        for i=1:n
            x(i) = range_min + (i-1) * delta;
            r2 = x(i);
        
            %compute switch scaling factor
            if(r2 <= switchon2)
                s = 1;
                ds = 0;
            end
        
            if(r2 > switchon2 && r2 < cutoff2)
                s = (cutoff2 - r2) * (cutoff2 - r2) * (cutoff2 + 2*r2 - 3 * switchon2) * denom;
                ds = 12 * (cutoff2 -r2) * (switchon2 - r2) * denom;
            end
        
            if(r2 >= cutoff2)
                s = 0;
                ds = 0;
            end
        
            inv_r = 1/(sqrt(r2));
            inv_r2 = 1/r2;
            inv_r6 = inv_r2 * inv_r2 * inv_r2;
            inv_r8 = inv_r6 * inv_r2;
            inv_r12 = inv_r6 * inv_r6;
            inv_r14 = inv_r12 * inv_r2;
        
            inv_r6_sw  = inv_r6 * s;
            inv_r12_sw = inv_r12 * s;      
        
            inv_r8_sw  = inv_r6  * (ds -  6*s*inv_r2);
            inv_r14_sw = inv_r12 * (ds - 12*s*inv_r2);
                 
            inv_r1 = cc * erfc(ewald_coef*sqrt(r2)) * inv_r;
            inv_r3 = -1 * inv_r2 * (cc * grad_coef * exp(-1*r2*ewald_coef2) + inv_r1);  
           
            fprintf(fileID1, '%d\t:\t%tx;\n', count, s);
            fprintf(fileID2, '%d\t:\t%tx;\n', count, ds);
            fprintf(fileID3, '%d\t:\t%tx;\n', count, inv_r6_sw);
            fprintf(fileID4, '%d\t:\t%tx;\n', count, inv_r12_sw);
            fprintf(fileID5, '%d\t:\t%tx;\n', count, inv_r8_sw);
            fprintf(fileID6, '%d\t:\t%tx;\n', count, inv_r14_sw);
            fprintf(fileID7, '%d\t:\t%tx;\n', count, inv_r1);
            fprintf(fileID8, '%d\t:\t%tx;\n', count, inv_r3);
        
            %fprintf(fileID1, '%5.20f\n', s);
            %fprintf(fileID2, '%5.20f\n', ds);
            %fprintf(fileID3, '%5.20f\n', inv_r6_sw);
            %fprintf(fileID4, '%5.20f\n', inv_r12_sw);
            %fprintf(fileID5, '%5.20f\n', inv_r8_sw);
            %fprintf(fileID6, '%5.20f\n', inv_r14_sw);
            %fprintf(fileID7, '%5.20f\n', inv_r1);
            %fprintf(fileID8, '%5.20f\n', inv_r3);
        
            count = count + 1;
        end
    
        range_min = range_min * 2;
        range_max = range_max * 2;
    end
    
    %p = polyfit(x,y,m);
    %f = polyval(p,x);
    %ncoef=length(p);
    
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);
    fclose(fileID5);
    fclose(fileID6);
    fclose(fileID7);
    fclose(fileID8);
    
    %x_val = 8+8/2048
    %f_val=polyval(p,x_val)
    %plot(x,y,'r.')
    %plot(x,y,'r.',x,f,'b.')
    