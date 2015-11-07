function col = t2c(temp,nc,cmap,cax)
%  col = t2c(temp,nc,cmap,cax)
%
% Converts a field to colors with 'nc' shades, using colormap 'cmap' and
% color axis 'cax'. Uses linear interpolation between contours.

% cmap can take the following values:
% - 'temp' (blue to red, with white contour straddling zero)
% - 'precip' (brown to green, with white contour straddling zero)
% - 'jet' (matlab's default heat map)
% - 'hsv' (hue-separated values)

% invoke colormap
if strmatch(cmap,'temp')
    cm = colormap(cejulien2(nc));
elseif strmatch(cmap,'precip')
    cm = colormap(ceprecipJEG(nc));
elseif strmatch(cmap,'jet')
    cm = colormap(jet(nc));
elseif strmatch(cmap,'hsv')
    cm = colormap(hsv(nc));
elseif strmatch(cmap,'hot')
    cm = colormap(flipud(hot(nc)));
elseif strmatch(cmap,'pmk')
    cm = colormap(pmkmp(nc));
elseif strmatch(cmap,'red')
    myColormap = [...
103,0,31
178,24,43
214,96,77
244,165,130
253,219,199
247,247,247]./255;
cm = flipud(myColormap); ncm = length(cm);
cx = linspace(1,ncm,nc);
cm = interp1(1:ncm,cm,cx);
    
end


tlims = linspace(cax(1),cax(2),nc);
n = length(temp);
col = zeros(n,3);
for k = 1:nc-1
    tidx = find(temp >= tlims(k) & temp < tlims(k+1));
    nidx = length(tidx);
    col(tidx,:) = repmat(cm(k,:),[nidx 1]);
end