@ ------------------------------- SolveDSS ----------------------------------

     Alfred Maußner
     28 January 2008
 
     Last Change: 03 April 2008

     This file contains a routine SolveVI used to compare to  Burkhard's
     value function iteration routine in RCh7_mont.g

     Contains routines that implement disreate state space methods
     for the solution of deterministic and stochastic DGE models.

     - MarkovAR
     - SolveVIS
     - _rhs_Bellman
     - _VI_valuefunction


------------------------------------------------------------------------------ @

@------------------------------------ MarkovAR ---------------------------------

   Alfred Maussner
   2004

   Purpose:  Approximate AR(1)-Process by Markov chain (Algorithm 9.2.1 of Heer and Maußner, 2005)
  
   Usage: {z,p}=MarkovAR(size,m,rho,sigma)
  
   Input: size:  scalar, the mulitple of the unconditional 
                 standard deviation of the AR(1) process
                 used to define the grid size
  
             m:  integer scalar, the number of grid points
   
           rho:  scalar, the autoregressive parameter of the process
  
         sigma:  scalar, the standard deviation of the innovations
  
   Output:   z:  m x 1 vector the grid approximating the process
  
             p:  m x m matrix of transition probabilities
----------------------------------------------------------------------------------- @
   
   
proc(2)=MarkovAR(size,m,rho,sigma);

  local p, zt, i, j, sigmaz, zbar;

  sigmaz=sqrt(sigma^2/(1-rho^2));
  zbar=size*sigmaz;  
  zt=seqa(-zbar,2*zbar/(m-1),m);

  p=zeros(m,m);
  i=1;
  Do until i>m;

   p[i,1]   =cdfn((zt[1]-rho*zt[i]+(zt[2]-zt[1])/2)/sigma); 
   
   j=2;
   do until j>(m-1);
      p[i,j]=cdfn((zt[j]-rho*zt[i]+(zt[j]-zt[j-1])/2)/sigma)-
             cdfn((zt[j]-rho*zt[i]-(zt[j]-zt[j-1])/2)/sigma);
           j=j+1;
   endo;
   p[i,m]=1-sumc(p[i,1:m-1]');
   i=i+1;
endo;
              
retp(zt,p);

endp;

@  --------------------------  SolveVIS -------------------------------------

   Purpose: Computes the policy function for a stochastic DGE model with
            one endogenous state variable x and one exogenous shock z.
                         
   Algorithm: The policy function is computed on a matrirx of of n by m points,
              where n is the number of points in the grid for x and m the
              number of points in the grid for z.

              The method employed is iteration over the value function.
              The user can choose interpolation between the grid points of x,
              either linear or cubic. For this to happen die gobal _VI_IP
              must be set either to 1 (for linear interpolation) or to
              2 (for cubic interpolation). The defaul is 0, i.e., no
              interpolation is used. In this case, the solution matrix xz
              is a matrix of intergers (i,j), where the (i,j) element
              of this matrix points to the index of the x that gives
              the optimal next-period value of the endogenous state variable.

              The  algorithm makes use of the concavity of the value function and the montonicity
              of the policy function. Specifically, we use a binary search algorithm to locate the maximum on the rhs of
              the Bellman equation.


      Usage:  {v1,xz}=SolveVIS(beta,xvec,zvec,pmat,v0);

  Input : beta:= the discount factor from the agent's life-time utility function.

          xvec:= n by 1 column vector with elements in ascending order that
                 represent valid elements of the variable x.

          zvec:= m by 1 column vector that represent the elements of the Markov chain
                 that represents the stochastic variable z

          pmat:= m by m matrix, the transition matrix of the Markov chain.

          v0  := n by m matrix with the elements of the inital value function.
                  
  Output:  v1 := n by m matrix, the value function at the solution

           xz  := n by m matrix,  if i=1, ..., n is the
                 index of the a point from xvec and j=1, 2, ..., m, is the index of point
                 in zvec, xy[i,j]=jstar is the index of xvec that gives the optimal
                 next-period value of the endogenous state variable x.

  Remarks:  The one-period return function rf(z,x,x') must be programmed in a procedure
            lhs=rf(z,x,x') with name rf and the specified ordering of the variables. lhs
            returns the one-period return, if the current state of the system is given by (x,z)
            and if a the next-period state x' is chosen.

-------------------------------------------------------------------------------------------------- @

/* The next lines set default values for the global variables used by the program */

declare matrix _VI_IP?=0;      @ =1 linear interpolation, =2, cubic interpolation, =0 without interpolation @
declare matrix _VI_MPI?=0;     @ =1 modified policy iteration, =0 without modified policy iteration @
declare matrix _VI_MPI_K?=30;  @ # number of iterations if _VI_MPI=1 @
declare matrix _VI_xvec?=0;  @ stores the x values used in interpolation @
declare matrix _VI_ymat?=0;  @ stores the y values uses in interpolation @
declare matrix _VI_zvec?=0;  @ stores the z values used in interpolation @
declare matrix _VI_pmat?=0;  @ stores the transition matrix used in interpolation @
declare matrix _VI_beta?=0;  @ stores information uses for linear interpolation @
declare matrix _VI_der?=0;   @ stores the derivative information used in cubic interpolation @ 
declare matrix _VI_xex?=0;   @ stores information to write the value function as a function of x[i] alone @
declare matrix _VI_zex?=0;   @ stores information to write the value function as a function of x[i] alone @

declare matrix _VI_eps?=0.01;  @ stopping criterium @
declare matrix _VI_nc?=50;     @ stop if number of consecutive iterations with unchanged indices of policy function exceeds this number @
declare matrix _VI_Max?=1000;  @ maximal number of intertions @
declare matrix _VI_BS?=1;      @ binaary search instead of sequential search @
declare matrix _VI_clear?=1;   @ clear screen for printing @
  
/* Here beings the code for Solve_VI_S */

proc(2)=SolveVIS(beta,xvec,zvec,pmat,v0);
                   
   local i, j, l, v1, v2, v3, h1, h2, t, nx, nz, dv, eps1, di, nc, w, js, jmin, jmax, jl, ju,
         ax, cx, bx, umat;

   external matrix _VI_MPI, _VI_MPI_K, _VI_xvec, _VI_ymat, _VI_xex, _VI_zex, _VI_eps, _VI_nc, _VI_der;
   external proc rf;

/* Step 1: Initialize */
eps1=_VI_eps*(1-beta);      @ convergence criteria @

_VI_beta=beta;

nx=rows(xvec);   @ number of grid points in xvec @
nz=rows(zvec);   @ number of grid points in zvec @

if _VI_IP==0; h1=ones(nx,nz); endif;            @ intial policy function @
h2=zeros(nx,nz); @ new policy function @ 
w =zeros(3,1);

v1=v0;              @ old policy function @
v2=zeros(nx,nz);    @ new policy function @
dv=1;
nc=0;

if _VI_IP/=0; _VI_xvec=xvec;_VI_zvec=zvec;_VI_ymat=v0; _VI_pmat=pmat; _VI_beta=beta; endif;
if _VI_IP==2; _VI_der=zeros(nx,nz); for j (1,nz,1); _VI_der[.,j]=CSpline(xvec,v1[.,j],1,0|0); endfor; endif;
if _VI_MPI==1; umat=zeros(nx,nz); v3=umat; endif;

/* Step 2: Iterate over the value function */
if _VI_clear; DosWinOpen("Value Function Iteration",0|0|15|1|1); cls; endif;
t=1;
do until (t>_VI_Max) or (dv<eps1) or (nc>_VI_nc); @ begin loop over value function @

    for j (1,nz,1);  @ begin loop over zvec @
        if _VI_IP/=0; _VI_zex=j; endif;
        js=1;
        for i (1,nx,1);  @ begin loop over xvec @
            if _VI_IP/=0; _VI_xex=xvec[i]; endif; 
            if _VI_BS;    
                jmin=js;                
                jmax=nx;
                do while (jmax-jmin)>2;       @ the next lines implement the binary search algorithm @
                    jl=floor((jmin+jmax)/2);
                    ju=jl+1;
                    w[1]=rf(zvec[j],xvec[i],xvec[jl])+beta*(pmat[j,.]*(v1[jl,.]'));
                    w[2]=rf(zvec[j],xvec[i],xvec[ju])+beta*(pmat[j,.]*(v1[ju,.]'));
                if w[2]>w[1]; jmin=jl; else; jmax=ju; endif;
                endo;
                w[1]=rf(zvec[j],xvec[i],xvec[jmin])+beta*(pmat[j,.]*(v1[jmin,.]'));
                if jmax>jmin;    w[2]=rf(zvec[j],xvec[i],xvec[jmin+1])+beta*(pmat[j,.]*(v1[jmin+1,.]')); else; w[2]=w[1]; endif;
                w[3]=rf(zvec[j],xvec[i],xvec[jmax])+beta*(pmat[j,.]*(v1[jmax,.]'));
                js=maxindc(w);
                if _VI_IP==0; v2[i,j]=w[js]; endif;               
                js=jmin+js-1;
               if _VI_MPI==1; umat[i,j]=rf(zvec[j],xvec[i],xvec[js]); endif;
            else;    
                jmin=js;
                w[1]=rf(zvec[j],xvec[i],xvec[jmin])+beta*(pmat[j,.]*(v1[jmin,.]'));         
                for jl (jmin+1,nx,1);                
                     w[2]=rf(zvec[j],xvec[i],xvec[jl])+beta*(pmat[j,.]*(v1[jl,.]'));
                    if w[2]<=w[1]; js=jl-1; if _VI_IP==0; v2[i,j]=w[1]; endif; if _VI_MPI==1; umat[i,j]=rf(zvec[j],xvec[i],xvec[js]); endif; break; else; w[1]=w[2]; endif;                
                endfor; 
            endif;
    
            /* The next lines implement linear interpolation between grid points */
            if _VI_IP/=0;
                if js==1;  /* boundary optimum, ax=bx=a[1]  */
                    ax=xvec[1];
                    bx=ax+eps1*(xvec[2]-xvec[1]);                      
                    cx=xvec[2];
                    if _rhs_bellman(j,xvec[i],bx)<_rhs_bellman(j,xvec[i],ax);
                        h2[i,j]=xvec[1];
                    else;                      
                        h2[i,j]=GSS(&_VI_valuefunction,xvec[1],xvec[2]);                                          
                    endif;
                elseif js==nx;   /* boundary optimum, bx=cx=a[n] */
                    ax=xvec[nx-1];
                    cx=xvec[nx];
                    bx=cx-eps1*(xvec[nx]-xvec[nx-1]);                    
                    if _rhs_bellman(j,xvec[i],bx)<_rhs_bellman(j,xvec[i],cx);
                        h2[i,j]=xvec[nx];
                    else;
                        h2[i,j]=GSS(&_VI_valuefunction,xvec[nx-1],xvec[nx]);                                            
                    endif;
                else;
                    h2[i,j]=GSS(&_VI_valuefunction,xvec[js-1],xvec[js+1]);                             
                endif;
    
                v2[i,j]=_rhs_bellman(j,xvec[i],h2[i,j]);    
            else;
                h2[i,j]=js;
            endif;
            
        endfor; @ end loop over xvec @
    endfor;    @ end loop over zvec @
    
    if _VI_IP==0;
        /* modified policy iteration */
        if _VI_MPI==1;            
            for l (1,_VI_MPI_K,1);                
                for i (1,nx,1);
                    for j (1,nz,1);
                        v3[i,j]=umat[i,j]+beta*(pmat[j,.]*(v2[h2[i,j],.]'));
                    endfor;
                endfor;
                v2=v3;
            endfor;
        endif;
        /* compute stopping criterium 2 */        
        di=sumc(sumc(h2 ./= h1));
        if di>=1; nc=0; else; nc=nc+1; endif; 
        h1=h2;
    endif;
    dv=maxc(maxc(abs(v2-v1))); 

    locate 5,5;
    ?"Iteration #= " ftos(t,"*.*lf",6,0);
    locate 6,5;
    ?"Largest element in v1-v0= " dv;
    if _VI_IP==0;
        locate 7,5; ?"# of indices that have changed: " di;
        locate 8,5; ?"# of consecutive iterations with constant policy function=" nc;
    endif;
    v1=v2;  
    if _VI_IP/=0; _VI_ymat=v1; endif;
    if _VI_IP==2; for j (1,nz,1); _VI_der[.,j]=CSpline(xvec,v1[.,j],1,0|0); endfor; endif;
    t=t+1;
endo;
if t>_VI_Max;
   locate 10,5;
   ?"Maximum number of iterations exceeded. Change _VI_Tmax!";
   locate 11,5;
   ?"The computed solution may be inaccurate.Press any key...";wait;
endif;
if _VI_IP==0;
    if minc(minc(h1))==1;  locate 12,5; ?"Policy function hits lower bound of grid"; endif;
    if maxc(maxc(h1))==nx; locate 13,5; ?"Policy function hits upper bound of grid"; endif;
else;
    if minc(minc(h2))==xvec[1];  locate 12,5; ?"Policy function hits lower bound of grid"; endif;
    if maxc(maxc(h2))==xvec[nx]; locate 12,5; ?"Policy function hits upper bound of grid"; endif;
endif;

retp(v2,h2);

endp;


@ ----------------------------------- _rhs_bellman ----------------------------------------------

  2 Feburary 2008
  Alfred Maussner

  Purpose: returns the rhs of the bellman equation rf(z0,x0,x1)+beta*pmat[z0,.]'v1[x1,.]'

  usage: vij=_rhs_bellman(z0,x0,x1);

  Input: z0 := scalar integer, the index of the current shock in the vector _VI_zvec;
         x0 := scalar, the current value of the endogenous state variable
         x1 := scalar, the next-period vaulue of the endogenous state variable

  Output: vij := ths value of the rhs of the Bellman equation at (z0,x0,x1)


  Remarks:   vij is found from interpolation (linear or cubic, depending on
             the global variable _VI_IP)
------------------------------------------------------------------------------------------------- @


proc(1)=_rhs_bellman(z0,x0,x1);

    local rhs, j, m;

    m=cols(_VI_pmat);
    rhs=rf(_VI_zvec[z0],x0,x1);
    if _VI_IP==1;
        for j (1,m,1);
            rhs=rhs+_VI_beta*_VI_pmat[z0,j]*LIP(_VI_xvec,_VI_ymat[.,j],x1);
        endfor;
    endif;
    if _VI_IP==2;
        for j (1,m,1);
            rhs=rhs+_VI_beta*_VI_pmat[z0,j]*Splint(_VI_xvec,_VI_ymat[.,j],_VI_der[.,j],1,x1);
        endfor;
    endif;

    retp(rhs);

endp;

@ ------------------------------------- _VI_valuefunction -----------------------------------------

   Alfred Maussner
   2 February 2008


   Purpose: The value function as a function of x1 alone. z0 and x0 are given in global variables.
            This is uses to find x' via golden section search.

   Usage:   vi=_VI_valuefunction(x1)

------------------------------------------------------------------------------------------------- @

proc(1)=_VI_valuefunction(x1);

    retp(_rhs_bellman(_VI_zex,_VI_xex,x1));

endp;

@ ---------------------------------------------- Ch4_Tools.src ---------------------------------------------------

    Alfred Mauï¿½ner
    
    Last Change: 08 September 2014
    
    Contains: Gauss procedures used by programs explained in Chapter 4
              of Heer/Mauï¿½ner, Dynamic General Equilibrium Modeling, 2nd. Ed.,
              Springer 2008.
              
   Content:
   
    GSS
    MyDate
    LIP
    BLIP
    CSpline
    Splint
    GH_Int4
    Macheps
    GraphSettings
    
------------------------------------------------------------------------------------------- @

struct mstr {string m;}; // structure use by MyDate
  struct mstr month;

@ ---------------------------------- GSS --------------------------------

   Alfred Maussner

   usage: x1=GSS(&f,xl,xu);

   Input: &f : pointer to the function f(x)
          xl : scalar, lower bound of the interval in which the maximizer
               lies
          xu : scalar, upper bound of the interval in which the maximizer
               lies.

   Output: x1, the approximate maximizer

------------------------------------------------------------------------- @

proc(1)=GSS(&f,xl,xu);

   local tol, p, q, a, b, c, d, fb, fc,i, f:proc;

  /* Compute the parameters of the problem */
   tol=sqrt(MachEps);

   p=(sqrt(5)-1)/2;
   q=1-p;

  /* Compute the initial interval [a,d] and the points b and c
  ** that divide it
 */
   a=xl; d=xu; b=p*a+q*d; c=q*a+p*d;

 /* Compute the function value at b and c */
   fb=f(b); fc=f(c);

 /* Iterate so that the inverval gets small and smaller */
   do while abs(d-a)>tol*maxc(1|(abs(a)+abs(c)));
      if fb<fc; /* choose [b,d] as next interval */
         a=b;
         b=c;
        fb=fc;
         c=p*b+q*d;
        fc=f(c);
      else; /* choose [a,c] as next interval */
         d=c;
         c=b;
        fc=fb;
         b=p*c + q*a;
         fb=f(b);
      endif;
   endo;

   if fb>fc; retp(b); else; retp(c); endif;

endp;


@ ------------------------------------- MyDate -------------------------------

   Purpose: Returns the current date and time from
            Gauss intrinsic commands and prints it to an open file.
          
            Day Month Year
            Hours:Minutes:Seconds

          

   Usage: MyDate;

--------------------------------------------------------------------------- @

proc(0)=MyDate;

  local t,d,str1,str2;
 
  month=reshape(month,12,1);

  month[1].m="January";
  month[2].m="Feburary";
  month[3].m="March";
  month[4].m="April";
  month[5].m="May";
  month[6].m="June";
  month[7].m="July";
  month[8].m="August";
  month[9].m="September";
  month[10].m="October";
  month[11].m="November";
  month[12].m="December";

  d=date;
  t=time;
  str1=ftos(d[3],"%0*.*lf",2,0) $+ " " $+ month[d[2]].m $+" " $+ ftos(d[1],"%0*.*lf",2,0);
  str2=ftos(t[1],"%*.*lf",2,0)  $+ ":" $+ ftos(t[2],"%0*.*lf",2,0) $+":" $+ ftos(t[3],"%0*.*lf",2,0);
  ?str1;
  ?str2;
retp;

endp;

@ --------------------------------- LIP --------------------------------------------

    25 January 2008
    Alfred Maussner

    Purpose: given a function y=f(x) tabulated in xvec and yvec and a point x0,
             return the function value y0=f(x0) obtained from linear interpolations
             between x1<x<x2.

    Usage: y0=LIP(xvec,yvec,x0);

   Input:  xvec: n by 1 vector, the tabulated values of the independent variable x
           yvec: n by 1 vector, the tabulated values of the dependent variable y
             x0: m by 1 vector, the values of the independent variable for which y0 is
                 to be computed.

   Output:   y0: m by 1 vector, see above.

------------------------------------------------------------------------------------ @

proc(1)=LIP(xvec,yvec,x0);

    local n, m, j, k, y0;
    
    n=rows(xvec);
    m=rows(x0);
    y0=zeros(m,1);

    for k (1,m,1);
        if (x0[k]<xvec[1]) or (x0[k]>xvec[n]); ?"Input out of grid. Procedure will return a missing value. Press any key"; wait; retp(miss(1,1)); endif;
        if x0[k]==xvec[1];
            y0[k]=yvec[1]; 
        elseif x0[k]==xvec[n];
            y0[k]=yvec[n];
        else;
            j=sumc(xvec.<=x0[k]); @ this determines the lower bracket for x0 @
        y0[k]=yvec[j]+((yvec[j+1]-yvec[j])/(xvec[j+1]-xvec[j]))*(x0[k]-xvec[j]);
        endif;
    endfor;
    retp(y0);
    
endp;

@ --------------------------------------- BLIP ------------------------------------------------

   2 February 2008
   Alfred Maussner
 
   Purpose: Bilinear Interpolation as given by Formulas (3.6.1) through (3.6.5)
            in Press et al. (1992), S. 116f.

   Usage: z=BLIP(xvec,yvec,zmat,x,y)
 
  Input: xvec = n by 1, vector, the grid of variable x
         yvec = m by 1, vector, the grid of variable y
         zmat = n by m matrix, with tabulated function values at
                      z(i,j)=f(x(i),y(j), i=1, ..., n, j=1, ...,m
            x = scalar,  the x coordinate
            y = scalar,  the y coordiante
 
  Output: z, the interpolated value of f(x,y)

  Remarks:  the elements in xvec and yvec must satisfy xvec[i]<xvec[i+1] for all i
            and similar for y
       
------------------------------------------------------------------------------------------------ @

proc(1)=BLIP(xvec,yvec,zmat,x,y);

    local n, m, i, j, z, t, u;

n=rows(xvec); m=rows(yvec);
z=zeros(4,1);

@ first, locate the square that surrounds (x,y) @
if (x<xvec[1]) or (x>xvec[n]); ?"x outside of grid! Program stops. Press any key";x;wait; retp(miss(1,1)); endif;
if (y<yvec[1]) or (y>yvec[m]); ?"y outside of grid! Program stops Press any key"; y;wait; retp(miss(1,1)); endif;

i=sumc(xvec.<=x);
j=sumc(yvec.<=y);

if (i==n) and (j==m);
    retp(zmat[n,m]);
elseif (i==n) and (j<m);
    u=(y-yvec[j])/(yvec[j+1]-yvec[j]);
    retp((1.-u)*zmat[n,j]+u*zmat[n,j+1]);
elseif (i<n) and (j==m);
    t=(x-xvec[i])/(xvec[i+1]-xvec[i]);
    retp(t*zmat[i+1,m]+(1.-t)*zmat[i,m]);
else;
    t=(x-xvec[i])/(xvec[i+1]-xvec[i]);
    u=(y-yvec[j])/(yvec[j+1]-yvec[j]);
    retp((1.-t)*(1.-u)*zmat[i,j]+t*(1.-u)*zmat[i+1,j]+t*u*zmat[i+1,j+1]+(1.-t)*u*zmat[i,j+1]);         
endif;
endp;


@ ------------------------------------ CSpline --------------------------------------------


   25 January 2008
   Alfred Maussner

   Purpose:   Computes the second derivatives for the CUBIC SPLINE APPROXIMATION of
              a function y=f(x). You must first run CSpline and then CSplint
              to evaluate the function.

   Useage:  y2=CSpline(x,y,cmethod,yp);


   Input:      x := n by 1 vector that stores the tabulated values of the independent variable
               y := n by 1 vector that stores the tabulated values of the dependent variable y=f(x)
          cmethod:= 1, 2, or 3:
                    if =1:  natural cubic spline, 
                    if =2:  secant hermite spline
                    if =3:  first derivaties at x[1] and x[n] are specified as yp[1] and yp[2]
               yp:= 2 by 1 vector, as just defined

   Output:     y2:= n by 1 vector, the second derivatives at x[i], i=1,n

   Remarks: The code is based on the Fortran subroutine spline in Press et al.(1992), p. 109.

------------------------------------------------------------------------------------------------- @

proc(1)=CSpline(x,y,cmethod,yp);

    local n, i, k, p, qn, sig, un, u, y2;

    /* Initializing */
    n=rows(x);
    u=zeros(n,1);
    y2=zeros(n,1);

    if cmethod==1;  @ natural cubic spline, i.e., y2[1]=y2[n]=0 @

            y2[1]=0.0;
             u[1]=0.0;
               qn=0.0;
               un=0.0;

    endif;

    if cmethod==2;  @ secant hermite spline, i.e., set the first derivative f' at x[1] and x[2] equal to the secant @

            yp[1]=(y[2]-y[1])/(x[2]-x[1]);
            yp[2]=(y[n]-y[n-1])/(x[n]-x[n-1]);
            y2[1]=-0.5;
            u[1]=(3./(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp[1]);
            qn=0.5;
            un=(3./(x[n]-x[n-1]))*(yp[2]-(y[n]-y[n-1])/(x[n]-x[n-1]));

    endif;

    if cmethod==3;  @ set second derivative to the values specified in yp @

            y2[1]=-0.5;
            u[1]=(3./(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp[1]);
            qn=0.5;
            un=(3./(x[n]-x[n-1]))*(yp[2]-(y[n]-y[n-1])/(x[n]-x[n-1]));

    endif;

    for i (2,n-1,1);

        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2;
        y2[i]=(sig-1.0)/p;
         u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;

    endfor;

    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   
    for k (n-1,1,-1);

        y2[k]=y2[k]*y2[k+1]+u[k];

    endfor;

    retp(y2);

endp;


@ --------------------------------------------------- CSplint ----------------------------------------


   25 January 2009
   Alfred Maussner
    
   Purpose: Computes the CUBIC SPLINE APPROXIMATION of a function y=f(x) being
            tabulated in the vectors x and y. y2 must be computed from x and y 
            via a first call to CSpline.

   Useage: call y0=Splint(x,y,y2,gmethod,x0);

   Input:       x := n by 1 vector, the tabulated values of the independent variable
                y := n by 1 vector, the tabulated values of the dependent variable
               y2 := n by 1 vector, the second derivative of f at x (computed from CSpline)
           gmethod:= 1 or 2
                     if =1: equally spaced nodes x[1]<x[2]<...<x[n]                     
                     if =2: bisection to locate x0[i]
               x0 := m by 1 vector, points at which the function is to be approximated

 Output:        y0:= m by 1 vector, the approximate values of f(x0);

 Remarks: Code is based on the Fortran routine splint from Press et al. (1992), p.110.

-------------------------------------------------------------------------------------------------- @

proc(1)=Splint(x,y,y2,gmethod,x0);

      local k,klo, khi, i, a, b, h, n, m, y0;

      n=rows(x);
      m=rows(x0);
     y0=zeros(m,1);
      
      for i (1,m,1);
         if (x0[i]<x[1]) or (x0[i]>x[n]); ?"x0 out of grid! The procedure stops and returns a missing value. Press any key";wait; retp(miss(1,1)); endif;
         if (x0[i]==x[1]);
             y0[i]=y[1];
         elseif (x0[i]==x[n]);
             y0[i]=y[n];
         else;
            if gmethod==1;
               klo=sumc(x.<=x0[i]);
               khi=klo+1;
            endif;
            if gmethod==2;  @ bisection @
                klo=1;
                khi=n;
                do while ((khi-klo)>.1);
                    k=(khi+klo)/2;
                    if (x[k]>x0[i]); khi=k;  else;  klo=k;  endif;
                endo;
            endif;
        
            h=x[khi]-x[klo];
            a=(x[khi]-x0[i])/h;
            b=(x0[i]-x[klo])/h;
           y0[i]=a*y[klo]+b*y[khi]+((a^3-a)*y2[klo]+(b^3-b)*y2[khi])*(h^2)/6;

         endif;

     endfor;

     retp(y0);

endp;


/*  GH_INT4: Gauss-Hermite integration over a one-dimensional 
**            space using four points
**
**   usage: x=GH_INT4(&f,sigma)
**
**  Input:      f, pointer to the function y=f(z), which is to be integrated
**          sigma, scalar, the standard deviation of the normal distribution
**
**      Output: x, value of the integral of y over [-inf,+inf]
**
**  Remark: The integration nodes and weights (are taken from Judd (1998), p. 262)
*/ 


proc(1)=GH_INT4(&f,sigma);

    local i, sum, s, temp, ghw, ghx, f: proc;

	ghx={-1.6506801230,-0.5246476232,0.5246476232,1.6506801230};
	ghw={0.08131283544,0.8049140900,0.8049140900,0.08131283544};

	sum=0.0;	
	s=sqrt(2.0)*sigma;
	for i (1,4,1);
			temp=f(s*ghx[i]);
			if scalmiss(temp); ?"Could not evaluate function! Press any key to continue";wait;retp(0); endif;
			sum=sum+temp*ghw[i];
	endfor;

    retp(sum/sqrt(pi));
	
endp;

@ ------------------------- MachEps ------------------------------

  Purpose: Computes the machine epsilon
  Usage: e=MachEps;

  Output: the smallest number e so that (1+e)>1 is true

------------------------------------------------------------------- @

proc(1)=MachEps;

    local eps, i;
        eps=1;
        do while (1+eps)>1;
            eps=eps/2;
        endo;
    retp(2*eps);
endp;

/* GraphSettings: Changes Gauss' default initialization of graphic routines
**
** Usage: GraphSettings;
*/
proc(0)=GraphSettings;

 external matrix  _pltype, _ptitlht, _pnumht,
         _paxht, _pmcolor, _plwidth, _pcsel, _pdate;
 external proc GraphSet;
 GraphSet;
   _pdate="";        @ do not plot date               @
  _pltype=6;         @ solid lines only               @
 _plwidth=4;         @ line width always 4            @
   _paxht=0.15;      @ size of axis labels            @        
  _pnumht=0.15;      @ size of axis numbering         @
 _ptitlht=0.20;      @ size of titles                 @
 _pmcolor={0,        @ color of axis: black           @
           0,        @ color of axis' numbers         @
           0,        @ color of x-axis label          @
           0,        @ color of y-axis label          @
           0,        @ color of z-axis label          @
           0,        @ color of title                 @
           0,        @ color of boxes                 @
           0,        @ color of date and author       @
           15};      @ background color:        white @
   _pcsel={0,        @ color of first line: black     @
           9,        @ color of second line: light blue @
           10,       @ color of third line:  light green @
           12,       @ color of fourth line: light red      @
           13,       @ color of fifth line:  light magenta @
           11,       @ color of sixth line: light cyan @
           6};       @ color of seventh line: brown @
   fonts("Simplex, Simgrma");
          
retp;
endp;           

