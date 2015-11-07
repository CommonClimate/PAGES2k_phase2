function [cs,h]=m_contourpm_white(long,lat,varargin);
%  M_CONTOURPM_WHITE Draws WHITE contour lines on a map
%    M_CONTOURPM(LONG,LAT,DATA,...) draw contours on a map. Behavior
%    is the same as for CONTOURPM except that LONG and LAT vectors or
%    matrices must be specified.
%
%    [CS,H]=M_CONTOURPM_WHITE(...) returns a contour matrix C and a vector
%    H of handles to LINE or PATCH objects for use by CLABEL.
%
%    See also CONTOUR
 %  
%   Julien Emile-Geay, Aug 11 2006. 
 % adapted from   m_contourpm by :
% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if min(size(long))==1 & min(size(lat))==1,
 [long,lat]=meshgrid(long,lat);
end;

[X,Y]=m_ll2xy(long,lat,'clip','on');
[cs,h]=contourpm_white(X,Y,varargin{:});
