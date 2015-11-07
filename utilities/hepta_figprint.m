function hepta_figprint(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function hepta_figprint(filename,res)
%
%   prints current figure to 'filename.eps' with resolution res
%    if no string is provided, abracadabra.eps is created.
%
%    res = resolution (dpi). Defalt = 800.
%
%   History : 
%        created by Hepta Technologies, GaTech, April 2007. 
%           Adapted from lprps by E. DiLorenzo
%	   modified by Hepta Technologies at NCAR, July 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
  filename = varargin{1};
  res= varargin{2}
elseif nargin == 1
  filename = varargin{1};
  res=800;
elseif nargin==0
  filename ='abracadabra.eps';
  res=800;
end

in=findstr(filename,'.eps');
if isempty(in)
	file_eps=[filename,'.eps'];
else
	file_eps=filename;
end	

%orient(gcf,'landscape','PaperPositionMode','auto');	
set(gcf,'PaperPositionMode','auto');
% in case .eps wasn't part of the name

% print to eps
disp(['EPS color file ...  ',file_eps]);
res_string=['-r' num2str(res)];
print('-depsc2', '-cmyk','-painters', '-loose', res_string , file_eps)
%



