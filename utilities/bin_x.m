function [Xb, tb] =  bin_x(t,X,bins,edgestring);
% function [Xb, tb ] =  bin_x(t,X,bins,edgestring);
%  Bin matrix X along axis t according to bins "bins"

if length(t) ~= size(X,1)
    error('age and value vectors must be the same length')
end

binstep = bins(2)-bins(1);
inclstart=0;
if nargin>3
    if strcmp(edgestring,'start')
        inclstart=1;
    end
end

for i = 1:length(bins)-1
    if inclstart
        q = find(t >= bins(i) & t<bins(i+1));
    else
        q = find(t > bins(i) & t<=bins(i+1));
    end
    Xb(i,:) = [nanmean(X(q,:),1)];
    tb(i) =  nanmean(t(q));
end

