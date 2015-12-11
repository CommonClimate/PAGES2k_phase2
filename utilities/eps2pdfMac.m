function eps2pdfMac(filename)
%eps to pdf
%uses system call to perl. 

if ~ismac
    warning('This has not been tested (and probably won''t work) on non-Mac operating systems; failure may be imminent')
end
%check for perl

if system(['unset LD_LIBRARY_PATH; export PATH="/Library/Frameworks/Python.framework/Versions/3.4/bin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/texbin"; perl -v']) ~=0
error('perl is either not installed on your system, or not in the PATH')
end

perlPdfPath=[pwd '/utilities/epstopdf.pl'];

system(['unset LD_LIBRARY_PATH; export PATH="/Library/Frameworks/Python.framework/Versions/3.4/bin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/texbin"; perl ' perlPdfPath ' ' filename]);
