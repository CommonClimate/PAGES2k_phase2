function [y]=ykskewness(data)
quantile25=quantile(data,0.25);
quantile50=quantile(data,0.5);
quantile75=quantile(data,0.75);
y=(quantile25-2*quantile50+quantile75)/iqr(data);
end