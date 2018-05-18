% Example:
% Compute coefficients for target function f=x^-7 on range from 2^0 to 2^8 for a 2nd or 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 or 256 intervals, and output results to a.txt.
%
% Command line:
% syms x
% f=x^-7
% poly_interpolation_coef_real(3,16,2^-4,2^-3,'a.txt')

function  LJ_poly_interpolation_coef(m,n,min,max,filename1, filename2)
% m is the order of interpolation. i.e, m=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if(nargin<6)
		error('LJ_poly_interpolation_coef(m,n,min,ma,filename1,filename2)');
	end
	if(min<0)
		error('min must be gt 0.');
	end

	fileID1 = fopen(filename1, 'wt');
    fileID2 = fopen(filename2, 'wt');
    
    cutoff = 12;
    switchon = 10;
    
    cutoff2 = cutoff * cutoff;
    switchon2 = switchon * switchon;
    inv_denom = (cutoff2 - switchon2)^3;
    denom = 1/inv_denom;
    count = 2048;
    
    x = zeros(n,1);
    y = zeros(n,1);
    
    ca = min;
    cb = min + (max - min)/n;
    
    delta = (max-min)/n;
    
    for i=1:n
        x(i) = min + (i-1) * delta;
        r2 = x(i);
        
        if(r2 <= 100)
            s = 1;
            ds = 0;
        end
        
        if(r2 > 100 && r2 <= 144)
            s = (cutoff2 - r2) * (cutoff2 - r2) * (cutoff2 + 2*r2 - 3 * switchon2) * denom;
            ds = 12 * (cutoff2 - r2) * (switchon2 - r2) * denom;
        end
        
        if(r2 > 144)
            s = 0;
            ds = 0;
        end
            
        y(i) = s;
        %fprintf(fileID, '%5.6f\n',s);
        fprintf(fileID1, '%d\t:\t%tx;\n', count, s);
        fprintf(fileID2, '%d\t:\t%tx;\n', count, ds);
        count = count + 1;
    end
    
    p = polyfit(x,y,m);
    f = polyval(p,x);
    ncoef=length(p);
    
    
    fclose(fileID1);
    fclose(fileID2);
    
    %x_val = 8+8/2048
    %f_val=polyval(p,x_val)
    plot(x,y,'r.')
    %plot(x,y,'r.',x,f,'b.')
    