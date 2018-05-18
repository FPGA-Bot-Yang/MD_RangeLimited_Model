% Example:
% Compute coefficients for target function f=x^-7 on range from 2^-4 to 2^7 for a 3rd order approximate polynomial with logirithmic szie intervals,
% every section is of 128 intervals, and output results to a.txt.
%
% Command line:
% syms x
% f=x^-7
% poly_interpolation_coef_real(3,16,2^-4,2^-3,'a.txt')

function  poly_interpolation_function(m,n,precision,min,max)
% m is the order of interpolation. i.e, m=1 produces ax+b
% the results are from lower order to higher order, i.e coef(0,0) is the coefficient of constant term for first bin.
	if(nargin<5)
		error('poly_interpolation_function(m,n,precision,min,max)');
	end
	if(min<0)
		error('min must be gt 0.');
	end

	fileID_0 = fopen('b0.txt', 'wt');
    fileID_1 = fopen('b1.txt', 'wt');
    fileID_2 = fopen('b2.txt', 'wt');
    fileID_3 = fopen('b3.txt', 'wt');
    
    count = 0;
    
    ewald_coef = 0.257952;
    ewald_coef2 = ewald_coef * ewald_coef;
    cc = 332.0636;
    grad_coef = 0.291067;
    
    range_min = min;
    range_max = 2*min;
    
    while(range_min < max)
    
        step = (range_max-range_min)/n;
        ca = range_min;
        %cb = range_min + step;
    
        for j=1:n
            x = zeros(precision,1);
            y = zeros(precision,1);
            u = zeros(precision,1);
            
            % for L-J Potential/Force
            r14 = zeros(precision,1);
            r8 = zeros(precision,1);
            
            delta = step/precision;
    
            for i=1:precision
                x(i) = ca + (i-1) * delta;
                k = x(i);
                
                inv_k = 1/k;
                inv_k2 = inv_k * inv_k;
                inv_k4 = inv_k2 * inv_k2;
                inv_k8 = inv_k4 * inv_k4;
                inv_k6 = inv_k2 * inv_k4;
                  
                %y(i) = -1 * cc * (1/k) *( erfc(ewald_coef*sqrt(k)) * (1/sqrt(k)) + grad_coef * exp(-1*k*ewald_coef2) );
                u(i) = cc * erfc(ewald_coef*sqrt(k)) * (1/sqrt(k));
                y(i) = -1 * (1/k) * ( cc * grad_coef * exp(-1*k*ewald_coef2) + u(i) );
               
                r8(i) = inv_k8;
                r14(i) = inv_k8 * inv_k6;
            end
    
            p = polyfit(x,y,m);
            q = polyfit(x,u,m);
            q = polyfit(x,r14,m);
            
            %f = polyval(p,x);
            ncoef=length(p);
        
            fprintf(fileID_0, '%d\t:\t',count);
            fprintf(fileID_1, '%d\t:\t',count);
            fprintf(fileID_2, '%d\t:\t',count);
            fprintf(fileID_3, '%d\t:\t',count);
        
            %fprintf(fileID_0,'%5.15f;',p(4));
            %fprintf(fileID_1,'%5.15f;',p(3));
            %fprintf(fileID_2,'%5.15f;',p(2));
            %fprintf(fileID_3,'%5.15f;',p(1));
            
            %fprintf(fileID_0,'%tx;',p(4));
            %fprintf(fileID_1,'%tx;',p(3));
            %fprintf(fileID_2,'%tx;',p(2));
            %fprintf(fileID_3,'%tx;',p(1));
            
            fprintf(fileID_0,'%tx;',q(4));
            fprintf(fileID_1,'%tx;',q(3));
            fprintf(fileID_2,'%tx;',q(2));
            fprintf(fileID_3,'%tx;',q(1));
            
            %fprintf(fileID_0,'%5.15f;',q(4));
            %fprintf(fileID_1,'%5.15f;',q(3));
            %fprintf(fileID_2,'%5.15f;',q(2));
            %fprintf(fileID_3,'%5.15f;',q(1));
            
            fprintf(fileID_0,'\n');
            fprintf(fileID_1,'\n');
            fprintf(fileID_2,'\n');
            fprintf(fileID_3,'\n');
        
            count = count + 1;
            ca = ca + step;
        end
        range_min = range_min * 2;
        range_max = range_max * 2;
    end
    
    fclose(fileID_0);
    fclose(fileID_1);
    fclose(fileID_2);
    fclose(fileID_3);
   
    