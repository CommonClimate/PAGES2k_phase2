function [istart,istop,segments] = find_breaks(t,x)
% FUNCTION [istart,istop,segments] = find_breaks(t,x)
%  identifies breaks in a timeseries interrupted with missing values
% (encoded by NaNs).  
%  Input: dependent variable x, independent variable t;
%  Output: istart, index of first non-NaN value;
%          istop, index of last non-NaN value;
%          segments, index sets of non-NaN values
%   JEG, USC, June 2012
% =======================================================	
nz=~isnan(x); dnz=diff(nz);

istart = find(nz==1,1,'first'); %ilast = find(nz==1,1,'last');
inext =  find(dnz==1)+1;
if (~isempty(inext) & sum(inext > istart)>0)
  jumps = find(dnz==1)+1;
  istart = unique([istart ; jumps]);
end
istop = unique([find(dnz==-1); length(nz)]);
nseg  = length(istart);

for i=1:nseg;
   segments{i}=[istart(i):istop(i)];
end

istart = min(istart);
istop = max(istop);
return
