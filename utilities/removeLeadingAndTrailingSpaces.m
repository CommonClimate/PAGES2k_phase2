function sout=removeLeadingAndTrailingSpaces(sin)


if iscell(sin)
    sout=sin;
    for i=1:length(sin)
        
        if ischar(sin{i}) %ignore if not a str
            sin2=sin{i};
            si=1;
            ei=length(sin2);
            while 1
                if  si==(ei+1)
                    break
                else
                    if ~strcmp(' ',sin2(si))
                        break
                    else
                        si=si+1;
                    end
                end
            end
            
            ei=length(sin2);
            while 1
                if ei==0
                    break
                else
                    if ~strcmp(' ',sin2(ei)) | ei==0
                        break
                    else
                        ei=ei-1;
                    end
                end
            end
        end
        if ei==0 | si==(length(sin)+1)
            sout{i}='NULL';
        else
            sout{i}=sin2(si:ei);
        end
    end
    
    
    
elseif ischar(sin) %ignore if not a str
    si=1;
    ei=length(sin);
    
    while 1
        if si==(ei+1)
            break
        else
            if ~strcmp(' ',sin(si))
                break
            else
                si=si+1;
            end
        end
    end
    
    while 1
        if ei==0
            break
        else
            if ~strcmp(' ',sin(ei))
                break
            else
                ei=ei-1;
            end
        end
    end
    
    if ei==0 | si==(length(sin)+1)
        sout='NULL';
    else
        sout=sin(si:ei);
    end
    
else
    sout=sin;
end
