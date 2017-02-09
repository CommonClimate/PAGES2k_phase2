function []=fancyplot_deco(title_str, xlabel_str, ylabel_str,varargin) 
% FUNCTION  []=fancyplot_deco(title_str, xlabel_str, ylabel_str,[FontSize]) 
%     Pretties up an otherwise plain-looking matlab plot
%   INPUTS:  - title (string)
%			    - xlabel (string)
%				 - ylabel (string)
%            - FontSize (default = 14), applies to title, other sizes are scaled down
%
%    Written 2008 by Julien Emile-Geay (GaTech)
% ================================================================
if nargin < 4 
	FontSize = 14;
else
	FontSize = varargin{1};
end

if nargin < 5 
	FontName = 'Helvetica';
else
	FontName = varargin{2};
end
	
hTitle=title(title_str);
hYLabel=ylabel(ylabel_str);
hXLabel=xlabel(xlabel_str);

set(gca,'FontName',FontName,'FontSize', round(FontSize*0.71));
set([hTitle, hXLabel, hYLabel], 'FontName' , FontName);
set([hXLabel, hYLabel],'FontSize', round(FontSize*0.86));
set(hTitle, 'FontSize', FontSize, 'FontWeight' , 'bold');
set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
