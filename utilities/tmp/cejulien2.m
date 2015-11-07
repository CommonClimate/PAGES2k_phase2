function [cmap] = cejulien2(nc);
% Function [cmap] = cejulien2(nc);
%
%    A very cute colormap with nc= 4*q+1 levels, q a small integer
%       (this is best for centered colorscales, with zero as white) 
 %  
%  INPUT  :   nc : small odd integer (13 or 17 recommended)
%   
%  OUTPUT :	cmap(nc,3)	- colormap RGB-value array
%     (c)    Julien Emile-Geay, LDEO, Aug 2006  
%            inspired from Gerd Krahman's ce(STUFF) series of colormaps  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set nc to default if not provide
if (nargin==0) 
	nc=13	
end		
   
% set 5 points
   % set 5 point map
C = [ 0.1324    0.1011    0.8799,
      0.3163    0.7411    1.0000,
      0.9941    1.0000    0.9882 ,
      1.0000    0.7783    0.1665 ,
      0.9358    0.0877    0.0419 ];
%     
nk=round((nc-1)/4);
cmap=zeros(nc,3);

for ki=1:4
   xo=[ki ki+1]; % old axis
   xn=linspace(ki,ki+1,nk+1); % new axis
   bounds=(ki-1)*nk+linspace(1,nk+1,nk+1);
   cmap(bounds,1)=interp1(xo,C(ki:ki+1,1),xn)';
   cmap(bounds,2)=interp1(xo,C(ki:ki+1,2),xn)';
   cmap(bounds,3)=interp1(xo,C(ki:ki+1,3),xn)';
end

