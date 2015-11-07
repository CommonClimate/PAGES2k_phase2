function [col,cm] = t2c_brewer(field,nc,scheme,cax, brighten_factor, flip)
%  col = t2c_brewer(temp,nc,cmap,cax)
%
% Converts a field to colors with 'nc' shades, using colormap 'cmap' and
% color axis 'cax'. Uses linear interpolation between contours.

% colormaps are generated using brewermap.m, which takes the arguments nc
% and scheme. 
%  See http://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--attractive-and-distinctive-colormaps
%  Hat tip to  Stephen Cobeldick for making such an awesome package
if nargin < 6
    flip = 0;
end
if flip
    cm = flipud(brewermap(nc,scheme));
else
    cm = brewermap(nc,scheme);    
end

if nargin >= 5
    disp('Brighten up!')
    cm = brighten(cm, brighten_factor);
end



tlims = linspace(cax(1),cax(2),nc);
n = length(field);
col = zeros(n,3);
for k = 1:nc-1
    tidx = find(field >= tlims(k) & field < tlims(k+1));
    nidx = length(tidx);
    col(tidx,:) = repmat(cm(k,:),[nidx 1]);
end


colormap(cm)
