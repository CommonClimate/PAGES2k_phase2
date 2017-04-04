function [cmin,cmax,cmp] = m_colorsig(ch,data,lon,lat,tag,mks, cax)
%
% function [cmin,cmax,cmp] = m_colorsig(ch,data,lon,lat,tag,mks,cax)
%
%   This function maps a colormap to station markers based on
%   correlation coefficient or another value. For use with m_map only. 
%   Map must be initialized using m_proj before calling m_colorsig.
%
%       input  :  ch      :       colormap handle
%                 data    :       station data vector
%                 lon     :       vector of longitudes
%                 lat     :       vector of latitudes
%                 tag     :       vector or number (if all same)
%                                 stationtype vector of tags:
%                                   1  .  point
%                                   2  o  circle
%                                   3  x  x-mark
%                                   4  v  triangle (down)
%                                   5  s  square
%                                   6  *  star
%                                   7  ^  triangle (up)
%                                   8  +  plus
%
%                 mks     :       [8] markersize (6 pt for filled shapes)
%
%                 cax     :       coloraxis (if different than your data)
%                                   [cmin cmax]
%
%       output :  cmin    :       colormin
%                 cmax    :       colormax
%                 cmp     :       colormap
%
%       uses   :  m_line.m
%
%   COLORMAP HANDLE is a single element to which the colormap has been assigned
%   STATION DATA is a vector that contains the meaning of the color
%   TAG is the same size as station_data, but is a tag to separate the types of
%       stations. starts from 1, can go up to 8. default is a vector of ones. 
%   MKS is the markersize for open shapes. Default is 8 for open shapes and 6
%       for filled shapes. The difference is always two.
%
%   version 0.1  J.Cherry july.2002  &  Felix Tubiana aug.2002
%
% function [cmin,cmax,cmp] = m_colorsig(ch,data,lon,lat,tag,mks,cax)
%
if nargin < 1, help m_colorsig, return, end
if nargin < 7, cax = []; end
if nargin < 6, mks = 8; end
if nargin < 5, tag = ones(length(data),1); end
if nargin < 4, error('not enough arguments'), end

if ~isempty(cax) & cax(2) <= cax(1)
   error('   caxis must be in the form [cmin cmax]')
end
if len(data) ~= len(lon) | len(data) ~= len(lat) 
   error('   Lon and Lat must be the same length as your data vector')
end
if len(tag) == 1
   tag = ones(len(data),1) * tag;
end
if len(tag) ~= len(data)
   error('   Tag vector must be the same length as your data vector')
end
if size(ch, 2) ~= 3 | max(ch(:)) > 1 | min(ch(:)) < 0
   error('   Please use a real colormap (n x 3)');
end
if max(tag) > 8 | min(tag) < 1
   error('   Please only use tag numbers 1 - 8.')
end

n = len(data);
m = min(data);

%% Normalize data to get a range between 0 and 1 %%
if isempty(cax)
   newdata = (data - m)/max(data - m);
   c = round(newdata.*(len(ch(1:end-1,1)))) + 1;
   % Use the normalized & rounded value index in the colormap

   cmp(1:n,:)  = ch(c(1:n,:),:);
   cmax = max(data);
   cmin = min(data);
   caxis([cmin cmax]);
else
   newdata = (data - cax(1))/(cax(2) - cax(1));
   newdata(newdata < 0) = 0;
   newdata(newdata > 1) = 1;
   c = round(newdata.*(len(ch(1:end-1,1)))) + 1;
   % Use the normalized & rounded value index in the colormap

   cmp(1:n,:)  = ch(c(1:n,:),:);
   cmax = cax(2);
   cmin = cax(1);
   caxis([cmin cmax]);
end

marker = [{'.'}; {'o'}; {'x'}; {'v'}; {'s'}; {'*'}; {'^'}; {'+'}];

% This alters the marker sizes so that they appear all the same size.

markersizes = ones(8,1)*mks - [-50 2 0 2 2 0 2 0 ]';
% Find all points with a given tag and plot them.
TAG = unique(tag);
for t = 1:len(TAG)
   pt = find(tag == TAG(t));
   for p = 1:len(pt)
      m_line(lon(pt(p)), lat(pt(p)), 'marker', char(marker(TAG(t))), ...
             'markersize', markersizes(TAG(t)), 'color', cmp(pt(p),:))
   end
end
