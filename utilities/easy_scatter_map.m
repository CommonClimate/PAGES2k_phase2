function easy_scatter_map(lonlat)


    m_proj('robinson');
    m_coast('patch',[0.5,0.5,0.5]);
    m_grid('linest','-','xticklabels',[],'yticklabels',[]);
    
    hold on;
    
    [X,Y] = m_ll2xy(lonlat(:,1), lonlat(:,2));
    
    scatter(X,Y,14,'filled');

end