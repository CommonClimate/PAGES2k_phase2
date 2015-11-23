function [ax,h1,h2] = yyplot(x,y1,y2,Xlab,Ylab,col)
% function [ax,h1,h2] = yyplot(t,X,Y,Xlab,Ylab,col)
%   Plot two variables, one as bars (X, h1) , the other as a curve (Y, h2),
%     as a function of time (t). 
%
% x : vector of length nx
% y1: vector of length nx
% y2: vector of length nx
% Ylab: cell array of strings {'Bar plot label',' Line plot label'}
% col : cell array {RGB row vector for y1, RGB row vector for y2}
% ========================================================================
style_l = {'FontName','Helvetica','FontSize',12,'FontWeight','bold'};

% plot solid line otherwise, and number of proxies
[ax,h1,h2]  = plotyy(x,y1,x,y2,@bar,@plot); hold on
set(h1,'edgecolor',col{1}); 
set(ax(1),'Ycolor',col{1},'YAxisLocation','right','TickDir','out','YMinorTick','on')
set(get(ax(1),'Ylabel'),'String',Ylab{1})
set(get(ax(1),'Ylabel'),style_l{:})
%set(h2,'visible','off'); %make invisible for now;
set(h2,'color',col{2},'linewidth',1);
set(ax(2),'Ycolor',col{2},'YAxisLocation','left')
set(get(ax(2),'Ylabel'),'String',Ylab{2});
set(get(ax(2),'Ylabel'),style_l{:})
linkaxes([ax(1),ax(2)],'x'); 
set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
xlabel(Xlab);
set(gca,style_l{:})
set(ax(2),'XTick',[])
ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
