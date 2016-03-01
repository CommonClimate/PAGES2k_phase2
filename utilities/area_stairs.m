function hhh = area_stairs(x,y)
%AREA_STAIRS  Filled area plot.
    %   AREA_STAIRS(X,Y) produces a stacked area plot suitable for showing the
    %   contributions of various components to a whole.
    %
    %   For vector X and Y, AREA_STAIRS(X,Y) is the same as STAIRS(X,Y) except that
    %   the area between 0 and Y is filled.  When Y is a matrix, AREA_STAIRS(X,Y)
    %   plots the columns of Y as filled areas.  For each X, the net
    %   result is the sum of corresponding values from the columns of Y.
    %
    %   See also AREA, PLOT, BAR.
    
    %   Published by: Hamed Amini (Hamed.Amini@hotmail.co.uk)
    %   Program by:   Oleg Komarov (oleg.komarov@hotmail.it)
    %   $Revision: 1.1 $  $Date: 20-Jul-2011$

if isvector(x) && ismatrix(y) && size(y,1)==length(x)
    x = reshape(x,length(x),1); 
    
else
    error('??? x and y dimentions should be consistent.')
end
    

x = [x';x'];

for i = 1:size(y,2)
    temp = [y(:,i)';y(:,i)'];
    temp = temp(:);
    Y(:,i) = temp;
end

hhh = area(x([2:end end]),Y);

end