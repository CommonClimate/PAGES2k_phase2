function xs = hepta_smooth(x,fc)
%  HEPTA_SMOOTH 
%   Lowpass filtering with AR padding on either end.
%   The padding length is taken as round(1/fc) where fc is the cutoff
%   frequency. Smoothing is achieved via a cubic spline as in [1] and [2].
%
%    Based on E.R.Cook's arpad.F (padding) and spline smoothing by 
%   [1] ER Cook, K Peters -  1981  Tree-Ring Bulletin, Vol. 41:
%       The Smoothing Spline: A New Approach To Standardizing
%       Forest Interior Tree-Ring Width Series For Dendroclimatic Studies
%   [2] Weinert, H. L., 2009. A fast compact algorithm for cubic spline
%   smoothing. Computational Statistics and Data Analysis 53, 932-940.
%
%   Note : Needs ARfit by Schneider & Neumaier, 2001
%     http://www.gps.caltech.edu/~tapio/arfit/ 
%
%      (Some Rights Reserved)  Hepta Technologies, May 2010      
%  ====================================================

% make sure x is a colum vector
[x,xm] = center(x(:)); n = length(x);
Tc = 1.0/fc;
np = ceil(Tc);
% Pad using AR predictions
xp = arpad(x,np,np);
% spline smoothing
amp = 0.5; % Amplitude reduction at period Tc. 
p      = 1/(((cos(2*pi*fc)+ 2)*(1-amp)/(12*amp*(cos(2*pi*fc)-1)^2)) + 1);
%lambda = p*Tc^3/(1-p); % Smoothing parameter in Weinert's formulation
[xs, score] = cspline(xp, 1, p); %Apply spline smoothing

xs = xs(np+1:np+n); % Restrict to original range
xs= xs(:)+xm; % force to column vetor

% ARPAD
% =============
function xp = arpad(x,nf,nb)
%  ARPAD 
%    Padding of a timeseries with ne points on either end.
%    Based on E.R.Cook's arpad.F
%   Use :
%
%     xp = arpad(x,nf,nb)
%
%    Inputs: - timseries x, length(n)
%            -  nf : number of forward pads
%            -  nb : number of back pads
%    Output: padded series xp 
%   v1.0 : Julien Emile-Geay, USC, May 18 2010
%   v1.1 : Julien Emile-Geay, USC, May 26 2010 
%      - adjusted maxmimum order of AR model
%

n = length(x);
% Put bounds on the AR order
p1 = max(nf,nb);
if (n-p1 <0)
   disp('Timeseries is too short for desired padding')
elseif p1 > round(n/5) % Heuristic Bound to ensure AR model stability
   p = round(n/5);
else
   p= p1;
end  
%
[w,A,C]=arfit(x,p,p); % AR fit 
%
xf = arpred(x(n-nf+1:n),A,nf); 
xu = arpred(flipud(x(1:nb)),A,nb); xb = flipud(xu);
%
xp  = [xb; x ; xf]; 

% AR Prediction
% =============
function [v]=arpred(x,A,np)
%ARPRED	Linear prediction of an AR process.	
%
%  v=ARPRED(x,A,C,np) simulates np time steps of the AR(p) process
%
%     v(k,:)' =   A1*x(k-1,:)' +...+ Ap*x(k-p,:)', 
%
%  where A=[A1 ... Ap] is the coefficient matrix, and x is a matrix.
%   x should be the last np rows of a matrix of interest.
%
% Inspired by arsim.m by Tapio Schneider (ARfit package)
%     www.gps.caltech.edu/~tapio/arfit/
%
%  Created May 6 2010
%    Julien Emile-Geay, USC
%     julieneg@usc.edu

  [n,m]       = size(x);               % dimension of state vectors 
  p       = size(A,2)/m;                % order of process

  if (p ~= round(p)) 
    error('Bad arguments.'); 
  end

  % Check whether specified model is stable
  A1 	  = [A; eye((p-1)*m) zeros((p-1)*m,m)];
  lambda  = eig(A1);
  if any(abs(lambda) > 1)
    warning('The specified AR model is unstable.')
  end
  
 
  % Get transpose of system matrix A (use transpose in simulation because 
  % we want to obtain the states as row vectors)
  AT      = A';

  mval = mean(x);
  z    = ones(p,1)*mval;
  
  % Initialize state vectors
  u      = [x; zeros(np,m)];
  
  % Simulate np observations. In order to be able to make use of
  % Matlab's vectorization capabilities, the cases p=1 and p>1 must be
  % treated separately.
  if p==1
    for k=2:np+1; 
      z(1,:) = u(k-1,:)*AT;
      u(k,:) = z;
    end
  else
    for k=p+1:np+p; 
      for j=1:p;
         z(j,:) = u(k-j,:)*AT((j-1)*m+1:j*m,:);
      end
      u(k,:) = sum(z);
    end
  end
  
  % return only the n simulated state vectors
  v = u(n+1:n+np,:); 



% Copy of H. Weinert's cpsline.m
% ==================================
function[x,score]=cspline(y,r,lam)
n=length(y); nc=ceil(n/2);
x=zeros(r*n-r+1,1); c=zeros(n,1); z=zeros(n,1);
w = y(1:n-2)-2*y(2:n-1)+y(3:n);
a0=6+lam*2/3; a1=4-lam/6;
    
    %Factor coefficient matrix, solve triangular systems, find trace
    
    e=zeros(1,n-2); f=zeros(1,n-2);
    d=a0; f(1)=1/d; c(2)=f(1)*w(1); mu=a1; e(1)=mu*f(1);
    d=a0-mu*e(1); f(2)=1/d; c(3)=f(2)*(w(2)+mu*c(2));
    mu=a1-e(1); e(2)=mu*f(2);
    for j=3:n-2
        d=a0-mu*e(j-1)-f(j-2); f(j)=1/d;
        c(j+1)=f(j)*(w(j)+mu*c(j)-c(j-1));
        mu=a1-e(j-1); e(j)=mu*f(j);
    end
    c(n-2)=c(n-2)+e(n-3)*c(n-1);
    for j=n-4:-1:1
        c(j+1)=c(j+1)+e(j)*c(j+2)-f(j)*c(j+3);
    end
    g2=f(n-2); tr1=g2; h=e(n-3)*g2; tr2=h;
    g1=f(n-3)+e(n-3)*h; tr1=tr1+g1; tr3=0;
    for k=n-4:-1:n-nc
        q=e(k)*h-f(k)*g2; tr3=tr3+q;
        h=e(k)*g1-f(k)*h; tr2=tr2+h; g2=g1;
        g1=f(k)*(1-q)+e(k)*h; tr1=tr1+g1;
    end
    q=e(n-nc-1)*h-f(n-nc-1)*g2; tr3=tr3+q;
    h=e(n-nc-1)*g1-f(n-nc-1)*h; tr2=tr2+h;
    tr1=6*(2*tr1-rem(n,2)*g1);
    tr2=-8*(2*tr2-(1+rem(n,2))*h);
    tr3=2*(2*tr3-rem(n,2)*q);
    tr=(tr1+tr2+tr3)/n;

%Compute GCV score

    z(1)=c(2); z(2)=c(3)-2*c(2);
    z(3:n-2) = c(2:n-3)-2*c(3:n-2)+c(4:n-1);
    z(n-1)=c(n-2)-2*c(n-1); z(n)=c(n-1);
    sq=(z'*z)/n; score=sq/tr^2;

%Compute estimates

if r > 1
    x(1:r:r*n-r+1) = y-z;
    for j = 1:r-1
        j1 = j/r;
        j2 = 1-j1;
        v = lam*j1*j2/6;
        j3 = v*(1+j1);
        j4 = v*(2-j1);
        for i = 1:n-1
            ir = i*r;
            x(ir-r+j+1) = j2*x(ir-r+1)+j1*x(ir+1)-j3*c(i+1)-j4*c(i);
        end
    end
else
    x = y-z;
end


