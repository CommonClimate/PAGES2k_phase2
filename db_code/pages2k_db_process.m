function  pages2k__db_process(vers,options)
% function pages2k__db_process_flat(vers)
%  INPUT:   vers (version of the database, string)
%           (optional) export. If true, exports figures to pdf. If not,
%           simply prints them to screen. [default =0]
%
% Tasks: Checks for duplicates, detrends corals, transforms to normality,
%  screens based on resolution & correlation

% load data and packages
year = [];
addpath(genpath('../utilities'))
fn = ['../data/PAGES2k_v' vers '_unpack.mat'];
load(fn)


%===========
%% transform to normality (this standardizes predictors as well)
% only used for reconstructions, not for QC plots
[nd,pd] = size(proxy_ann);
proxy_n = NaN([nd,pd]);


%% Visual check the effect of normalization, i.e. are we making up non-sense data?
warning off;
alpha = 0.05;
disp('Test for normality...')
for k = 1:pd
    X = proxy_ann(:,k); X = X(~isnan(X)); 
    if isempty(X)
        error([S(k).paleoData_pages2kID, ' (record',int2str(k),') is just NaNs'])
    end
    % test for normality
    [H_sp(k), pValue_sp(k)] = spiegel_test(X,alpha); % Spiegelhalter
    [H_li(k), pValue_li(k)] = lillietest(X,alpha);   % Lilliefors
    [H_jb(k), pValue_jb(k)] = jbtest(X,alpha);       % Jarque-Bera

    % test for symmetry
    [H_sy(k), pValue_sy(k)] = symmetry_test(X,alpha); % Lanzante [1996]
    %
    % different measures of skewness
    S1(k)  = (max(X)- median(X))/(median(X)- min(X)); % "B" statistic
    NPS(k) = (mean(X) - median(X))/std(X);  % non-parametric skewness [-1;1]
    YKS(k) = ykskewness(X);   % Yule-Kendall Skewness
end

disp('Gaussianize records...')
for j = 1:pd
    noNaN = ~isnan(proxy_ann(:,j));
    proxy_n(noNaN,j)  = gaussianize(proxy_ann(noNaN,j));
end
save(fn,'proxy_n','-append','-v6')

% identify as suspect those proxies that fail either the Spiegelhalter or
% the symmetry test
skewed  = find(H_sy); N_skew = numel(skewed);


if options.export
    %Plot all histograms and records
    lab{1} = 'Raw';
    lab{2} = 'ITS';

    fig(1), step = 3;
    set(gcf,'position',[224   234   1056   572]);
    disp('Plot gaussianized records...')
    for n=1:step:pd
        clf
        for k = n:n+step-1
            if k <=pd
                X = standardize(proxy_ann(:,k));
                nz = ~isnan(X);

                % title
                rec = ['#' num2str(k)];
                obj = archive{k};

                lat = p_lat(k);
                if (lat>=0)
                    lat_str=[num2str(round(lat)),'{}^{\circ}N'];
                else
                    lat_str=[num2str(-round(lat)),'{}^{\circ}S'];
                end

                lon = p_lon(k);
                if (lon>=0)
                    lon_str=[num2str(round(lon)),'{}^{\circ}E'];
                else
                    lon_str=[num2str(-round(lon)),'{}^{\circ}W'];
                end
                site_n = strrep(recordNames{k}, '_', '\_');
                %str = [rec,': ', obj,', ', site_n,' (',lon_str,',',lat_str,')'];
                str = [rec,': ', obj,', ', site_n];

                % plot timeseries
                subplot(3,3,k-n+1)  % change this
                plot(year(nz),X(nz),'color',Graph{p_code(k),1}), hold on;
                plot(year(nz),proxy_n(nz,k),'color',rgb('Silver'));
                fancyplot_deco(str,'Time','Transformed proxy',11);
                axis tight
                hl = legend(lab{:},'location','best'); set(hl,'box','off')
                % plot histogram
                subplot(3,3,k-n+4)
                histfit(X(nz));
                h = findobj(gca,'Type','patch');
                set(h,'FaceColor',Graph{p_code(k),1});
                h = findobj(gca,'Type','line');
                set(h,'Color',rgb('Black'));
                ylabel('PDF','FontName','Times','FontSize',10)
                xlabel('x','FontName','Times','FontSize',10)
                str=[rec,': ', obj,', ', site_n,', RAW'];

                title(str,'FontName','Times','FontSize',10)

                %  write out result of normality and symmetry test
                if H_sp(k)
                    nor = 'Non-normal';
                else
                    nor = 'Normal';
                end
                if H_sy(k)
                    sy = 'Skewed';
                else
                    sy = 'Symmetric';
                end

                text(.7,.8,['p = ', sprintf('%0.2e',pValue_sp(k)),',' nor],'units','normalized');
                text(.7,.7,['p = ', sprintf('%0.2e',pValue_sy(k)),',' sy ],'units','normalized');
                text(.7,.6,['YKS = ', sprintf('%+3.2f',YKS(k))],'units','normalized');
                text(.7,.5,['NPS = ', sprintf('%+3.2f',NPS(k))],'units','normalized');

                subplot(3,3,k-n+7)
                histfit(proxy_n(nz,k));
                h = findobj(gca,'Type','patch');
                set(h,'FaceColor',rgb('Silver'));
                h = findobj(gca,'Type','line');
                set(h,'Color',rgb('Black'))
                ylabel('PDF','FontName','Times','FontSize',10)
                xlabel('x','FontName','Times','FontSize',10)
                str=[rec,': ', obj,', ', site_n,', ITS'];
                title(str,'FontName','Times','FontSize',10)
            end
        end
        if n < pd
            filen=['../../figs/skewness/symtest0p01_', num2str(n),'_to_',num2str(n+step-1),'.pdf'];
        elseif n == pd
            filen=['../../figs/skewness/symtest0p01_', num2str(n),'.pdf'];
        end
        %if options.export
        export_fig(filen,'-cmyk','-r300','-nocrop');
        %end
    end
end

save(fn,'-append','NPS','YKS','H_sy','H_sp','H_li','H_jb','-v6')

% Q: Test again for normality? Those records that are still not normal after being Gaussianized will be booted out unceremoniously
% A:  everything is "Gaussian" after that (as far as we can tell...)

if  options.InterpSuperAnn
    superAnn = find(resMed>1);
    figDir = '../figs/annualize';
    proxy_a   = pages2k_db_annualizer(proxy_ann,S,superAnn,year,resMed,options.export,figDir);
    proxy_na  = pages2k_db_annualizer(proxy_n   ,S,superAnn,year,resMed);
    %proxy_da  = pages2k_db_annualizer(proxy_d   ,S,superAnn,year,resMed);
    %proxy_nda = pages2k_db_annualizer(proxy_nd  ,S,superAnn,year,resMed);
end

save(fn,'-append','proxy_a','proxy_na','-v6')
