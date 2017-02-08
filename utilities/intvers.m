function i=intvers;
% function INTVERS
% usage    i=intvers
% returns the leading integer of the running MATLAB version for tests 
% like 'if intversion>3; ...'

% Uwe Send, IfM Kiel, Mar 1993

vers=version;
i=eval(vers(1));
