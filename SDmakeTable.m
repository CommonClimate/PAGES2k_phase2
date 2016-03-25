%   make table for the SD paper

S = pages2k.S; nr = length(S);

dataSetName={S.dataSetName}';
latValue=({S.geo_meanLat})';   % NEEDS ROUNDING
longValue=({S.geo_meanLon})';
%
archiveType=pages2k.archive';
% chronological information
minYear=round(arrayfun(@(x) min(x.year),S)');
maxYear=round(arrayfun(@(x) max(x.year),S)');
minYear(find(minYear==0))=-1;
maxYear(find(maxYear==0))=-1;
minYear=num2cell(minYear);
maxYear=num2cell(maxYear);
yearUnits=({S.yearUnits})';

elevation=({S.geo_elevation})';
siteName=({S.geo_siteName})';
geo_pages2kRegion=({S.geo_pages2kRegion})';


proxyType=({S.paleoData_proxy})';
parameter=({S.paleoData_parameter})';
desc = {S.paleoData_description}';


climateInterpretation_climateVariable=({S.climateInterpretation_climateVariable})';
climateInterpretation_climateVariableDetail=({S.climateInterpretation_climateVariableDetail})';
climateInterpretation_seasonality=({S.climateInterpretation_seasonality})';
climateInterpretation_basis=({S.climateInterpretation_basis})';
climateInterpretationDirection=({S.climateInterpretation_interpDirection})';

% transfer to numeric format
sgn = lower({S.climateInterpretation_interpDirection});
sgn_vec = zeros(nr,1);
for r = 1:nr
    if strcmpi(sgn{r},'positive') || strcmpi(sgn{r},'p')  % KLUDGE
        sgn_vec(r) = +1;
    elseif strcmpi(sgn{r},'negative')
        sgn_vec(r) = -1;
    end
end


ID=({S.paleoData_TSid})';
authors=({S.pub1_author})';
investigators=({S.investigators})';
pubYear=({S.pub1_year})';
DOI=({S.pub1_DOI})';
QCnotes={S.paleoData_QCnotes}';

for ii=1:length(investigators)
    if size(investigators{ii,1},2)>0
        if isstr(investigators{ii,1})
            investigators2{ii,1}=investigators{ii,1};
        else
            investigators2{ii,1}=investigators{ii,1}(1);
        end
    end
end
noAuthor=find(cellfun(@length,authors)==0);
yesInvestigator=find(cellfun(@length,investigators2)>0);
repAuthor=intersect(noAuthor,yesInvestigator);
authors(repAuthor)=investigators2(repAuthor);
% extract first author
for ii=1:nr
    aName=authors{ii};
    if ~isempty(aName)
        if ~isstr(aName)
            aName=authors{ii}{1};
        end
        if isstr(aName)
            strend=min([strfind(aName,';') strfind(aName,',') strfind(aName,' ')]);
            if ~isempty(strend)
                firstAuthor{ii,1}=aName(1:(strend-1));
            else
                firstAuthor{ii,1}=aName;
            end
        end
    end
    %pubCite{ii} = [firstAuthor{ii,1} 'et al.,[' int2str(pubYear{ii}) ']']
end

for ii=1:length(siteName)
    sName=siteName{ii};
    if ~isempty(sName)
        if ~isstr(sName)
            sName=siteName{ii}{1};
        end
        sName(ismember(sName,' ,.:;!-')) = [];
        shortSite{ii,1}=sName(1:min(length(sName),12));
    end
end
for ii=1:length(siteName)
    year=pubYear{ii};
    if ~isempty(year)
        pYear{ii,1}=year(1:4);
    end
end


dataSetName={S.dataSetName}';

for ii=1:length(DOI)
    doi=cell2str(DOI{ii,1});
    doi(ismember(doi,'}{,;!')) = [];
    DOI{ii,1}=doi;
end
for ii=1:length(authors)
    as=cell2str(authors{ii,1});
    as(ismember(as,'}''{')) = [];
    authors{ii,1}=as;
end
for ii=1:length(siteName)
    sN=cell2str(siteName{ii,1});
    sN(ismember(sN,'}''{')) = [];
    siteName{ii,1}=sN;
end

%screenFDR={S.screenFDR}';
%signifHRneighbors=num2cell(arrayfun(@(x) sum(x.signifHRneighbors),S)');
%badS=arrayfun(@(x) isempty(x.signifHRneighbors),S)';
%signifHRneighbors(badS)={''};
medRes = [S.resMed];

medResolution=cellstr(num2str(round(medRes)'));
subann = find(medRes<1); ns = length(subann);
medResolution(subann) = {'<1'};

%
inComposite = cellstr(repmat('no',[nr 1]));
inComposite(idx_q) = {'yes'}; 
%%
%repeat for all regions
header={'ID' 'Name' 'Lat' 'Lon' 'Archive' 'Proxy' 'minYear' 'maxYear' 'Resolution' 'Include in composite' 'First Author' 'Publication Year' 'DOI' };
outcell=[ID dataSetName  latValue  longValue archiveType proxyType minYear maxYear medResolution inComposite firstAuthor pYear DOI];
cell2csv(['./data/SD-Table_v' vers '.csv'],[header ; outcell],',');
%
% notcell=~arrayfun(@(X) iscell(X.paleoData_values),TS);
% unSum=nan(length(TS),1);
% unSum(notcell)=arrayfun(@(X) nansum(X.paleoData_values),TS(notcell));
%
