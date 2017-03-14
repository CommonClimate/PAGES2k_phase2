function [] = saveMe(filename,varargin)

if ismac
    save(filename,varargin,'-v6')
else
    save(filename,varargin)
end