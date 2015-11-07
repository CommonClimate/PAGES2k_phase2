function H=fig(varargin)
% H=fig([name/number][,figureproperties])
%
% 
% sets these default props: 'NumberTitle','off','PaperOrientation','landscape'
% 
% Aslak Grinsted 2003

figs=get(0,'children');

if length(varargin)>0
    if ischar(varargin{1})
        sName=varargin{1};
        figs=figs(find(strcmp(get(figs,'type'),'figure')));
        figs=figs(find(strcmp(get(figs,'name'),sName)));
        if length(figs)>0
            H=figure(figs(1));
        else
            H=figure;
        end
        set(H,'name',sName);
    else
        sName='';
        number=varargin{1};
        H=figure(varargin{1});
        set(H,'name',['Fig ' num2str(H)]);
    end
else
    H=figure;
    set(H,'name',['Fig ' num2str(H)]);
end

if ~any(figs==H)
    set(H,'NumberTitle','off');
    orient(H,'landscape','PaperPositionMode','auto');
end
if length(varargin)>1
    set(H,varargin{2:end});
end
set(H,'color',[1 1 1])

if nargout==0
    clear H
end
