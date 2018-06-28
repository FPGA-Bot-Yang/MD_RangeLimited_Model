% Example:
% Compute coefficients for target function f=x^-7 on range from 2^-4 to 2^7 for a 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 intervals, and output results to a.txt.
%
% Command line:
% syms x
% f=x^-7
% poly_interpolation_coef_real(3,16,2^-4,2^-3,'a.txt')

function  poly_interpolation_coef(m,n,min,max,filename)
% m is the order of interpolation. i.e, m=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if(nargin<5)
		error('poly_interpolation_coef_real(m,n,min,ma,filename)');
	end
	if(min<0)
		error('min must be gt 0.');
	end

	fileID = fopen(filename, 'wt');
    
    ewald_coef = 0.257952;
    ewald_coef2 = ewald_coef * ewald_coef;
    cc = 332.0636;
    grad_coef = 0.291067;
    
    x = zeros(n,1);
    y = zeros(n,1);
    
    ca = min;
    cb = min + (max - min)/n;
    
    delta = (max-min)/n;
    
    for i=1:n
        x(i) = min + (i-1) * delta;
        k = x(i);
        y(i) = -1 * cc * (1/k) *(erfc(ewald_coef*sqrt(k)) * (1/sqrt(k)) + grad_coef * exp(-1*k*ewald_coef2));
    end
    
    p = polyfit(x,y,m)
    %f = polyval(p,x);
    ncoef=length(p);
    
    fclose(fileID);
    
    %x_val = 8+8/2048
    %f_val=polyval(p,x_val)
    %plot(x,f,'r.')
    %plot(x,y,'r.',x,f,'b.')
    