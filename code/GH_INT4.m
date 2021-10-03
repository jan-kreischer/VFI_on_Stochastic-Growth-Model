function integ = GH_INT4(f,sigma)
%  GH_INT4:    Gauss-Hermite integration over a one-dimensional space using four points
%  Input:      f, pointer to the function y=f(z), which is to be integrated
%              sigma, scalar, the standard deviation of the normal distribution
%  Output:     x, value of the integral of y over [-inf,+inf]
%
%  Remark:     The integration nodes and weights (are taken from Judd (1998), p. 262)

    ghx=[-1.6506801230;-0.5246476232;0.5246476232;1.6506801230];
	ghw=[0.08131283544;0.8049140900;0.8049140900;0.08131283544];
	sum=0.0;	
	s=sqrt(2.0)*sigma;
    for i=1:4
			temp=f(s*ghx(i));
            if ismissing(temp)
                disp("Could not evaluate function! Press any key to continue");
                pause;
                integ=missing;
                return
            end
			sum=sum+temp*ghw(i);
    end
    integ=sum/sqrt(pi);
end

