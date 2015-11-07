function [maximum, minimum] = maxmin(data_array)
% function [maximum, minimum] = maxmin(data_array)
%
%   Returns the max and min.
%
%   Displays maximum, minimum, mean, and Std. Dev. over all dimensions.
%

% Felix Tubiana
%    5.14.02
 
if nargin ~=1
   help maxmin
   return
end
A = reshape(data_array, numel(data_array), 1);
maximum = max(A);
minimum = min(A);
disp([' ']);
disp(['   Minimum Value: ',num2str(minimum),'   Maximum Value: ', ...
       num2str(maximum)]);
end
