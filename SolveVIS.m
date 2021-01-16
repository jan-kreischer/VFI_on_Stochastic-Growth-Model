function [v2,h2] = SolveVIS(beta,xvec,zvec,pmat,v0)

%Purpose: Computes the policy function for a stochastic DGE model with
%         one endogenous state variable x and one exogenous shock z.
%                          
%Algorithm: The policy function is computed on a matrirx of of n by m points,
%           where n is the number of points in the grid for x and m the
%           number of points in the grid for z.
% 
%           The method used: iteration over the value function.
%           There are three options for the interpolation between the grid
%           points of x: 
%           1. linear (_VI_IP = 1) 
%           2. cubic  (_VI_IP = 2)
%           3. no interpolation (_VI_IP = 0)
%
%           In this case, the solution matrix xz
%           is a matrix of intergers (i,j), where the (i,j) element
%           of this matrix points to the index of the x that gives
%           the optimal next-period value of the endogenous state variable.
% 
%           The  algorithm makes use of the concavity of the value function and the montonicity
%           of the policy function. Specifically, we use a binary search algorithm to locate the maximum on the rhs of
%           the Bellman equation.
% 
%Input : beta:= the discount factor from the agent's life-time utility function.
% 
%           xvec:= n by 1 column vector with elements in ascending order that
%                  represent valid elements of the variable x.
% 
%           zvec:= m by 1 column vector that represent the elements of the Markov chain
%                  that represents the stochastic variable z
% 
%           pmat:= m by m matrix, the transition matrix of the Markov chain.
% 
%           v0  := n by m matrix with the elements of the inital value function.
%                   
%Output:  v1 := n by m matrix, the value function at the solution
% 
%         xz  := n by m matrix,  if i=1, ..., n is the
%                index of the a point from xvec and j=1, 2, ..., m, is the index of point
%                in zvec, xy[i,j]=jstar is the index of xvec that gives the optimal
%                next-period value of the endogenous state variable x.
% 
%Remarks:  The one-period return function rf(z,x,x') must be programmed in a procedure
%             lhs=rf(z,x,x') with name rf and the specified ordering of the variables. lhs
%             returns the one-period return, if the current state of the system is given by (x,z)
%             and if a the next-period state x' is chosen.


    global VI_eps VI_beta VI_IP VI_xvec VI_zvec VI_ymat VI_pmat VI_xex VI_zex VI_Max VI_nc VI_BS delta
    eps1=VI_eps*(1-beta);      % convergence criteria 
    
    VI_beta=beta;
    
    nx=length(xvec);   % number of grid points in xvec 
    nz=length(zvec);   % number of grid points in zvec   
    
    h2=zeros(nx,nz);  % new policy function 
    w=zeros(3,1);
    
    v1=v0;              % old policy function 
    v2=zeros(nx,nz);    % new policy function 
    dv=1;
    nc=0;
    
    if VI_IP == 1
        VI_xvec = xvec;
        VI_zvec = zvec;
        VI_ymat = v0;
        VI_pmat = pmat;
        VI_beta = beta;
    end

%begin loop over value function      
t=1;
while (t<=VI_Max) && (dv>=eps1) && (nc<=VI_nc)
  
    for j = 1:nz    % begin loop over zvec
        if VI_IP==1 
            VI_zex=j;
        end
        js=1;
        for i=1:nx   % begin loop over xvec
            if VI_IP==1
                VI_xex=xvec(i); 
            end
            if VI_BS
                jmin=js;
                while xvec(jmin)<(1-delta)*xvec(i)
                    jmin=jmin+1;
                end                                
                jmax=nx;
                while (jmax-jmin)>2      % the binary search algorithm
                    jl=floor((jmin+jmax)/2);
                    ju=jl+1;
                    c1=rf(zvec(j),xvec(i),xvec(jl));
                    c2=rf(zvec(j),xvec(i),xvec(ju));
                    w(1)=c1+beta*(pmat(j,:)*(v1(jl,:)')); 
                    w(2)=c2+beta*(pmat(j,:)*(v1(ju,:)'));
                    if w(2)>w(1)
                        jmin=jl;
                    else
                        jmax=ju;
                    end
                end
                w(1)=rf(zvec(j),xvec(i),xvec(jmin))+beta*(pmat(j,:)*(v1(jmin,:)'));
                if jmax>jmin
                    w(2)=rf(zvec(j),xvec(i),xvec(jmin+1))+beta*(pmat(j,:)*(v1(jmin+1,:)'));
                else
                    w(2)=w(1);
                end
                w(3)=rf(zvec(j),xvec(i),xvec(jmax))+beta*(pmat(j,:)*(v1(jmax,:)'));
                [~,js]=max(w);
                if VI_IP==0
                    v2(i,j)=w(js);
                end              
                js=jmin+js-1;
                
            else 
                jmin=js;
                w(1)=rf(zvec(j),xvec(i),xvec(jmin))+beta*(pmat(j,:)*(v1(jmin,:)'));         
                for jl=jmin+1:nx                
                     w(2)=rf(zvec(j),xvec(i),xvec(jl))+beta*(pmat(j,:)*(v1(jl,:)'));
                    if w(2)<=w(1)
                        js=jl-1; 
                    else
                        w(1)=w(2); 
                    end
                end
            end   
            
            % Implementation of linear interpolation between grid points
            if VI_IP==1
                if js==1  % boundary optimum
                    ax=xvec(1);
                    bx=ax+eps1*(xvec(2)-xvec(1));
                    cx=xvec(2);
                    if rhs_bellman(j,xvec(i),bx)<rhs_bellman(j,xvec(i),ax)
                        h2(i,j)=xvec(1);
                    else                   
                        h2(i,j)=GSS(@VI_valuefunction,xvec(1),xvec(2));                                          
                    end
                elseif js==nx   % boundary optimum
                    ax=xvec(nx-1);
                    cx=xvec(nx);
                    bx=cx-eps1*(xvec(nx)-xvec(nx-1));
                    if bx>(1-delta)*xvec(i)
                        if rhs_bellman(j,xvec(i),bx)<rhs_bellman(j,xvec(i),cx)
                            h2(i,j)=xvec(nx);
                        else
                            if xvec(nx-1)>=(1-delta)*xvec(i)
                                h2(i,j)=GSS(@VI_valuefunction,xvec(nx-1),xvec(nx));
                            else
                                h2(i,j)=GSS(@VI_valuefunction,(1-delta)*xvec(i),xvec(nx));
                            end
                        end
                    else
                    h2(i,j)=xvec(nx);
                    end
                else
                    if xvec(js-1)>=(1-delta)*xvec(i)
                        h2(i,j)=GSS(@VI_valuefunction,xvec(js-1),xvec(js+1));
                    else
                        h2(i,j)=GSS(@VI_valuefunction,(1-delta)*xvec(i),xvec(js+1));
                    end
                end            
                v2(i,j)=rhs_bellman(j,xvec(i),h2(i,j));    
            else
                h2(i,j)=js;
            end
            
        end % end loop over xvec
  
    end     % end loop over zvec
    
%     if VI_IP==0
%         % compute stopping criterium 2       
%         di=sum(sum(h2~=h1));
%         if di>=1
%             nc=0;
%         else
%             nc=nc+1;
%         end
%         h1=h2;
%     end
    dv=max(max(abs(v2-v1))); 
    
    clc
    fprintf("Largest increment in V %f\n", dv)
    fprintf("Iteration #= %d\n", t);
    
%     if VI_IP==0
%         disp("# of indices that have changed: ");
%         disp(di)
%         disp("# of consecutive iterations with constant policy function=");
%         disp(nc);
%     end
    
    v1=v2;  
    if VI_IP==1
        VI_ymat=v1;
    end
    t=t+1;
end
if t>VI_Max
   disp("Maximum number of iterations exceeded. Change VI_Max!");
   disp("The computed solution may be inaccurate.Press any key...");
   pause;
end
if VI_IP==0
    if min(min(h1))==1
        disp("Policy function hits lower bound of grid");
    end
    if max(max(h1))==nx
        disp("Policy function hits upper bound of grid");
    end
else
    if min(min(h2))==xvec(1)
        disp("Policy function hits lower bound of grid");
    end
    if max(max(h2))==xvec(nx)
        disp("Policy function hits upper bound of grid");
    end
end
end
