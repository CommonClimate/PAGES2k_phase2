function h = m_line_fewer(long,lat,num,varargin);
% M_LINE_FWER Create a line on a map, with fewer markers
%    M_LINE(LONG,LAT,NUM) adds the line in vectors LONG and LAT to the
%    current axes. If LONG and LAT are matrices the same size, one
%    line per column is added.  NUm is the number of markers
%    
%
%    LINE returns a column vector of handles to LINE objects,
%    one handle per line. LINEs are children of AXES objects.
%
%    The LONG,LAT pair can be followed by
%    parameter/value pairs to specify additional properties of the lines.
%    These are standard 'line' properties.
%
%    See also LINE, M_LL2XY


% Julien Emile-Geay, USC, 02-Jan-2017
% hybrid of Richa Pawlocwicz's m_line.m (M-map package), with line_fewer_makers.m by Massimo Ciacci
% http://www.mathworks.com/matlabcentral/fileexchange/42560-line-fewer-markers

clp='on'; num = 5;

k=1;
while k<length(varargin),
    switch lower(varargin{k}(1:3)),
        case 'cli',
            clp=varargin{k+1};
            if isempty(findstr(clp,'on')),
                varargin{k+1}='off';
            else
                varargin{k+1}='on';
                %        varargin([k k+1])=[];
            end;
            k=k+2;
        otherwise
            k=k+2;
    end;
end;

[X,Y]=m_ll2xy(long,lat,'clipping',clp);

if nargout>0,
    [h1 h2 h3] = line_fewer_markers(X,Y,num,varargin{:});
else
    line_fewer_markers(X,Y,num,varargin{:});
end;


