function [h]=m_pcolor(long,lat,data,varargin);
%  M_PCOLOR Makes a pcolor image on a map.
%    M_PCOLOR(LONG,LAT,DATA,...) is a pseudocolor  plot of matrix DATA.
%    The values of the elements of DATA specify the color in each
%    cell of the plot. In the default shading mode, 'faceted',
%    each cell has a constant color and the last row and column of
%    DATA are not used. With shading('interp'), each cell has color
%    resulting from bilinear interpolation of the color at its 
%    four vertices and all elements of DATA are used. 
%    The smallest and largest elements of DATA are assigned the first and
%    last colors given in the color table; colors for the remainder of the 
%    elements in DATA are determined by table-lookup within the remainder of 
%    the color table.
%
%    See also M_CONTOUR, CONTOURF, PCOLOR

% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 19/02/98 - type - should have been 'clip','patch', rather than 'off'.
%  9/12/98 - handle all-NaN plots without letting contour crash.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)


global MAP_PROJECTION 

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if min(size(long))==1 & min(size(lat))==1,
 [long,lat]=meshgrid(long,lat);
end;


if abs(long(2,1)-long(1,1)) > abs(long(1,2)-long(1,1)) % x increases downward
     long=long';
     lat=lat';
     data=data';
end


dx=diff(long(1,1:2));
dy=diff(lat(1:2,1));
long=[long; long(end,:)]; long=[long long(:,end)+dx]-dx/2; 
lat=[lat; lat(end,:)+dy]; lat=[lat lat(:,end)]-dy/2; 

% duplicate last row and column so that pcolor will show them
data=[data; data(end,:)]; data=[data data(:,end)]; 


[X,Y]=m_ll2xy(long,lat,'clip','on');  %First find the points outside

i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;

                 % And then recompute positions without clipping. THis
                 % is necessary otherwise contouring fails (X/Y with NaN
                 % is a no-no. Note that this only clips properly down
                 % columns of long/lat - not across rows. In general this
                 % means patches may nto line up properly a right/left edges.
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','patch'); end;  

if any(~i(:)),
 [h]=pcolor2(X,Y,data,varargin{:});
 set(h,'tag','m_pcolor');
else
  h=[];
end;

if nargout==0,
 clear  h
end;


function h=pcolor2(x,y,z)
% pcolor2(x,y,z): pcolor front end that plots everything and assumes
% x and y are the mid-point location2 of the grid boxes.
%
% Example:
%   [x,y]=meshgrid(3:7,1:4); z=rand(size(x));
%   subplot(2,1,1), pcolor(x,y,z), shading flat
%   subplot(2,1,2), pcolor2(x,y,z)
%
% Ian Eisenman, 2006

% if nargin<3 % x and y not given
%     % remove singleton dimensions
%     x=squeeze(x);
%     z=x;
%     [x,y]=meshgrid(1:size(z,2),1:size(z,1));
% end
% 
% if min(size(x))==1 % x,y are vectors, not matrices
%     [x,y]=meshgrid(x,y);
% end
% 
% % remove singleton dimensions
% z=squeeze(z);
% 
% % make sure x is rows and y is columns
% if abs(x(2,1)-x(1,1)) > abs(x(1,2)-x(1,1)) % x increases downward
%     x=x';
%     y=y';
%     z=z';
% end
% 
% dx=diff(x(1,1:2));
% dy=diff(y(1:2,1));
% x=[x; x(end,:)]; x=[x x(:,end)+dx]-dx/2; 
% y=[y; y(end,:)+dy]; y=[y y(:,end)]-dy/2; 
% 
% % duplicate last row and column so that pcolor will show them
% z=[z; z(end,:)]; z=[z z(:,end)]; 

h0=pcolor(x,y,z);

% shading flat;

if(nargout > 0)
    h = h0;
end

