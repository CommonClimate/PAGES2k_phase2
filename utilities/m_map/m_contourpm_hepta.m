function [cs,h]=m_contourpm_hepta(long,lat,varargin);
%  M_CONTOURPM_HEPTA(LONG,LAT,DATA,CVEC,zlw,pmlw,cmat)  Draws  contour lines on a map
%    Behavior is the same as for CONTOURPM except that LONG and LAT vectors or
%    matrices must be specified.
%
%    [CS,H]=M_CONTOURPM_HEPTA(...) returns a contour matrix C and a vector
%    H of handles to LINE or PATCH objects for use by CLABEL.
%
%    See  CONTOURPM_HEPTA
%  
%   Hepta Technologies, Feb 2007. 
%   ---------------------------------------------
% adapted from  m_contourpm by Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%   requires ; M_Map package, Hepta Package.
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
[cs,h]=contourpm_hepta(X,Y,varargin{:});
