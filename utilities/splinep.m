function p = splinep(per,amp)
% splinep:  spline smoothing parameter from Cook's period of spline response
% CALL: p = splinep(per,amp);
%
% Meko 3-4-92
%
%************  IN 
%
% per (1 X 1)r - period in years at which amplitude of freq response should be amp
% amp (1 x 1)r - desired amplitude of freq response at period per years
%
%************ OUT 
%
% p (1 x 1)r spline smoothing parameter, following Cook and Peters (1981) 
%  corresponding to a spline with the desired amplitude of frequency response
%  at period per years.  This p is same as the "spline smoothing parameter" P
%  required as input to the MATLAB spline toolbox function csaps.m.
%___________________________________________________________________

% Convert period per in frequency f

f=1.0/per;

% Use eq. (2) of Cook and Peters to find p

p=6*amp*(cos(2*pi*f)-1)^2/((cos(2*pi*f)+2)*(1-amp));


% End of file
