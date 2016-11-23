% SDmakeTable(vers)  make excel table for the SD paper

nr = length(S);

archive = pages2k.archive;
dataSetName={S.dataSetName}';
latValue=({S.geo_meanLat})';   % NEEDS ROUNDING
longValue=({S.geo_meanLon})';
%
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
desc = {S.paleoData_description}';


climateInterpretation_variable=({S.climateInterpretation_variable})';
climateInterpretation_variableDetail=({S.climateInterpretation_variableDetail})';
climateInterpretation_seasonality=({S.climateInterpretation_seasonality})';
climateInterpretation_basis=({S.climateInterpretation_basis})';
climateInterpretationDirection=({S.climateInterpretation_interpDirection})';

% % transfer to numeric format
% sgn = cellfunlower({S.climateInterpretation_interpDirection});
% sgn_vec = zeros(nr,1);
% for r = 1:nr
%     if strcmpi(sgn{r},'positive') || strcmpi(sgn{r},'p')  % KLUDGE
%         sgn_vec(r) = +1;
%     elseif strcmpi(sgn{r},'negative')
%         sgn_vec(r) = -1;
%     end
% end


ID=({S.paleoData_pages2kID})';
[IDs isort] = sort(ID);

authors=({S.pub1_author})';
investigators=({S.investigators})';
pubYear=({S.pub1_pubYear})';

citeKey=({S.pub1_citeKey})';
noRef = cellfun(@isempty,citeKey); % locus of missing references
citeKey(noRef) = {'MISSING'};      % flag missing references


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
            strend=min([strfind(aName,';') strfind(aName,',')]);
            %strend=min([strfind(aName,';') strfind(aName,',') strfind(aName,' ')]);

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
        if ischar(year)
            pYear{ii,1}=year(1:4);
        else
            pYear{ii,1} = int2str(year);
        end
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
medRes = [pages2k.resMed];

medResolution=cellstr(num2str(round(medRes)));
subann = find(medRes<1); ns = length(subann);
medResolution(subann) = {'<1'};

%
inComposite = cellstr(repmat('no',[nr 1]));
inComposite(idx_q) = {'yes'};

% proxy-proxy relations
psig = zeros(nr,1);
signif_n = pages2k.signif_n;
for r = 1:nr
    if ~isempty(signif_n{r})
        psig(r) = sum(signif_n{r});
    end
end

%%  WRITE OUT FILES
% 1) csv table
header={'ID' 'Name' 'Lat' 'Lon' 'Archive' 'Proxy' 'minYear' 'maxYear' 'Resolution' 'First Author' 'Publication Year' 'CiteKey' 'inComposites'};
outcell=[ID(isort) dataSetName(isort)  latValue(isort)  longValue(isort) archive(isort)' proxyType(isort) minYear(isort) maxYear(isort) medResolution(isort) firstAuthor(isort) pYear(isort) citeKey(isort) inComposite(isort)];
cell2csv(['./SD-Table_v' vers '.csv'],[header ; outcell],',');
% 2) BibTeX citations
citeKey_c = unique(citeKey(~noRef)); nrefs = length(citeKey_c)

refStr = '\cite{'; % initialize cite command
for r = 1:nrefs-1
   refStr = [refStr,citeKey_c{r},',']; % fill in citekeys
end
refStr = [refStr, citeKey_c{nrefs}, '}']; % finalize  cite command
% save to latex file
%read_write_entire_textfile('PAGES2k_refs.tex', refStr)

%  Legacy debug hacks
% =====================
% notcell=~arrayfun(@(X) iscell(X.paleoData_values),TS);
% unSum=nan(length(TS),1);
% unSum(notcell)=arrayfun(@(X) nansum(X.paleoData_values),TS(notcell));

%find pages2k_phase1 paper

%phase1 = strmatch('ahmed2013continental',citeKey)
%dataSetName(phase1)

%
