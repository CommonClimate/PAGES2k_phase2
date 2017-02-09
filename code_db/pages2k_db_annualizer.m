function proxy_a = pages2k_db_annualizer(proxy,S,index,year,resMed,varargin);
%  PAGES2K_DB_ANNUALIZER : annualize super-annual records using a combination
%  of interpolation techniques.
%  Small, infrequent gaps: interpolation via singular spectrum analysis
%  Large gaps:    nearest-neighbor interpolation 
%   proxy_a = pages2k_db_annualizer(proxy,S,index,year,ResMed,[export=0],[figdir]);
%
%
%  INPUT:  - proxy: proxy matrix of size (nt x nr)
%          - S:  database structure (length nr)
%          - index : index of records that need to be annualized (length <= nr)
%          - year  : yearly time axis (length nt)
%          - ResMed; vector of median resolution (length nr)
%          - export: boolean (default = 0): if yes, export figures
%          - figDir: figure directory.
%
%  OUTPUT  - proxy_a : proxy, with annalized data in the 'index' columns
%
% History: Written by Julien Emile-Geay at USC, Aug 28 2015.

numvarargs = length(varargin);

if numvarargs == 0
    export = 0;
elseif numvarargs == 1
    export = varargin{1};
    figDir = './';
elseif numvarargs == 2
    export = varargin{1};
    figDir = varargin{2};
else
    disp('Optional arguments beyond 2 will be scorned')
end

nrecs = numel(find(index));
op.crit_name = 'pct_var'; op.P = 1.00; % SSA parameters

% initialize proxy matrix
proxy_a = proxy;
% then loop through records and overwrite columns corresponding to 'index'
for s  = 1:nrecs
    r = index(s);
    T = S(r); tc = year;
    disp(['Annualizing record #' int2str(s) '/' int2str(nrecs) ': ', T.dataSetName ]), clf
    
    Xc = proxy(:,r);
    noNaNx = ~isnan(Xc);
    tn = tc(noNaNx); Xn = Xc(noNaNx);
    ta = [min(tn):max(tn)]';  dtn = diff(tn);
    if resMed(r) <= 5
        F = griddedInterpolant(tn,Xn,'linear'); Xa = F(ta); % regrid first
        [spec,eig_vec,PC,RC,Xann,modes] = hepta_ssam(Xa,op); % SSA interpolation
        
    else   % otherwise, just assign to closest year, and remove potential duplicates
        F = griddedInterpolant(tn,Xn,'nearest'); Xann = F(ta);
    end
    
    if export
        fig('Interpolation'), clf
        %plot(ta,Xann,'color',Graph{p_code(r),1},'linewidth',2); hold on;
        plot(ta,Xann,'r-','linewidth',2); hold on;
        plot(tn,Xn,'marker','o','color',rgb('Gray'),'MarkerSize',8,'linestyle','none');
        %xlim(min(ta),max(ta)
        hold off; lab{2} = 'Raw'; lab{1} = 'Annualized'; legend(lab{:}), legend boxoff
        ylab = [T.paleoData_variableName ' (' removeLeadingAndTrailingSpaces(T.paleoData_units) ')'];
        site_n = strrep(T.dataSetName,'_','\_');
        ttl = [int2str(r),') ', T.archiveType,', ', site_n];
        fancyplot_deco(ttl,'Year (CE)',ylab,14);
        %pause
        filen=[figDir '/annualize/pages_2k_phase2_record_', sprintf('%03d',r) '.pdf'];
        export_fig(filen,'-r100','-cmyk','-painters','-nocrop')
    end
    
    % export to proxy matrix
    proxy_a(ismember(year,ta),r) = Xann;
end

end
