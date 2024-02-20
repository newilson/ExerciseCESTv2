function out = NWexerciseCESTsort(in)
%
%
% Guesses scan directories in exercise CEST protocol
%
% Use 'pre' and 'post' for pre and post exercise
% Use 'ref' or 'none' for reference scan


if nargin>0
    scandirs = dir(in);
else
    scandirs = dir;
end

fold = strsplit(in,'\');
fold = fold{end};

cestprecount = 0;
cestpostcount = 0;
out = []; 
for zz = 1:2 % loop through 2 times
    for ii=1:length(scandirs)
        if scandirs(ii).isdir
            lowname = lower(scandirs(ii).name);
            splitname1 = strsplit(lowname,'_');
            if strfind(fold,'RESEARCH_DJ') % CHOP
                lowname1 = lowname;
            else
                if length(splitname1)>3
                    lowname1 = strjoin(splitname1(end-3:end),'_'); % ignores sequence name part
                else
                    lowname1 = [];
                end
            end
            if strfind(lowname,'b0') % B0
                if strfind(lowname,'post')
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0post1','var')
                            B0post2 = ind;
                        else
                            B0post1 = ind;
                        end
                    elseif isequal(B0post1+1,B0post2)
                        if isequal(ind,B0post1)
                            out.B0magpost = scandirs(ii).name;
                        elseif isequal(ind,B0post2)
                            out.B0phpost = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    elseif isequal(B0post2+1,B0post1)
                        if isequal(ind,B0post2)
                            out.B0magpost = scandirs(ii).name;
                        elseif isequal(ind,B0post1)
                            out.B0phpost = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    else
                        disp([ii,zz])
                        warning('something unexpected')
                    end
                elseif strfind(lowname,'_pre')
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0pre1','var')
                            B0pre2 = ind;
                        else
                            B0pre1 = ind;
                        end
                    elseif isequal(B0pre1+1,B0pre2)
                        if isequal(ind,B0pre1)
                            out.B0magpre = scandirs(ii).name;
                        elseif isequal(ind,B0pre2)
                            out.B0phpre = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    elseif isequal(B0pre2+1,B0pre1)
                        if isequal(ind,B0pre2)
                            out.B0magpre = scandirs(ii).name;
                        elseif isequal(ind,B0pre1)
                            out.B0phpre = scandirs(ii).name;
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    else
                        disp([ii,zz])
                        warning('something unexpected')
                    end
                else
                    splitname = strsplit(lowname,'_');
                    if strcmp(lowname(1),'s') && ~isnan(str2double(lowname(2:5)))
                        ind = str2double(lowname(2:5));
                    else
                        ind = str2double(splitname(end));
                    end
                    if zz==1
                        if exist('B0unk1','var')
                            B0unk2 = ind;
                        else
                            B0unk1 = ind;
                        end
                    else
                        if isequal(B0unk1+1,B0unk2)
                            if ~isfield(out,'B0magpost') && isfield(out,'B0magpre')
                                if isequal(ind,B0unk1)
                                    out.B0magpost = scandirs(ii).name;
                                elseif isequal(ind,B0unk2)
                                    out.B0phpost = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            elseif ~isfield(out,'B0pre') && isfield(out,'B0post')
                                if isequal(ind,B0unk1)
                                    out.B0magpre = scandirs(ii).name;
                                elseif isequal(ind,B0unk2)
                                    out.B0phpre = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            end
                        elseif isequal(B0unk2+1,B0unk1)
                            if ~isfield(out,'B0magpost') && isfield(out,'B0magpre')
                                if isequal(ind,B0unk2)
                                    out.B0magpost = scandirs(ii).name;
                                elseif isequal(ind,B0unk1)
                                    out.B0phpost = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            elseif ~isfield(out,'B0pre') && isfield(out,'B0post')
                                if isequal(ind,B0unk2)
                                    out.B0magpre = scandirs(ii).name;
                                elseif isequal(ind,B0unk1)
                                    out.B0phpre = scandirs(ii).name;
                                else
                                    disp([ii,zz])
                                    warning('something unexpected')
                                end
                            end
                        else
                            disp([ii,zz])
                            warning('something unexpected')
                        end
                    end
                end
            elseif strfind(lowname1,'b1') % B1
                if strfind(lowname,'post')
                    out.B1post = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.B1pre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out,'B1post') && isfield(out,'B1pre')
                            out.B1post = scandirs(ii).name;
                        elseif ~isfield(out,'B1pre') && isfield(out,'B1post')
                            out.B1pre = scandirs(ii).name;
                        else
                            out.unknown.B1 = scandirs(ii).name;
                        end
                    end
                end
            elseif strfind(lowname1,'cest') % CEST images
                if strfind(lowname,'post')
                    out.CESTpost = scandirs(ii).name;
                    cestpostcount = cestpostcount+1/2;
                elseif strfind(lowname,'_pre')
                    out.CESTpre = scandirs(ii).name;
                    cestprecount = cestprecount+1/2;
                else
                    if zz==2
                        if ~isfield(out,'CESTpost') && isfield(out,'CESTpre')
                            out.CESTpost = scandirs(ii).name;
                            cestpostcount = cestpostcount+1;
                        elseif ~isfield(out,'CESTpre') && isfield(out,'CESTpost')
                            out.CESTpre = scandirs(ii).name;
                            cestprecount = cestprecount+1;
                        else
                            out.unknown.CEST = scandirs(ii).name;
                        end
                    end
                end
            elseif strfind(lowname1,'wassr') % WASSR images
                if strfind(lowname,'post')
                    out.WASSRpost = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.WASSRpre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out,'WASSRpost') && isfield(out,'WASSRpre')
                            out.WASSRpost = scandirs(ii).name;
                        elseif ~isfield(out,'WASSRpre') && isfield(out,'WASSRpost')
                            out.WASSRpre = scandirs(ii).name;
                        else
                            out.unknown.WASSR = scandirs(ii).name;
                        end
                    end
                end
            elseif ~(isempty(strfind(lowname1,'none')) && isempty(strfind(lowname1,'ref'))) % Reference scan
                if strfind(lowname,'post')
                    out.Refpost = scandirs(ii).name;
                elseif strfind(lowname,'_pre')
                    out.Refpre = scandirs(ii).name;
                else
                    if zz==2
                        if ~isfield(out,'Refpost') && isfield(out,'Refpre')
                            out.Refpost = scandirs(ii).name;
                        elseif ~isfield(out,'Refpre') && isfield(out,'Refpost')
                            out.Refpre = scandirs(ii).name;
                        else
                            out.unknown.Ref = scandirs(ii).name;
                        end
                    end
                end
            end
        end
    end
end

if ~isfield(out,'Refpost'), out.Refpost = ''; end
if ~isfield(out,'Refpre'), out.Refpre = ''; end
if ~isfield(out,'WASSRpost'), out.WASSRpost = ''; end
if ~isfield(out,'WASSRpre'), out.WASSRpre = ''; end
if ~isfield(out,'CESTpost'), out.CESTpost = ''; end
if ~isfield(out,'CESTpre'), out.CESTpre = ''; end
if ~isfield(out,'B1post'), out.B1post = ''; end
if ~isfield(out,'B1pre'), out.B1pre = ''; end
if ~isfield(out,'Refpost'), out.Refpost = []; end
if ~isfield(out,'B0magpost') || ~isfield(out,'B0phpost'), out.B0magpost = ''; out.B0phpost = ''; end
if ~isfield(out,'B0magpre') || ~isfield(out,'B0phpre'), out.B0magpre = ''; out.B0phpre = ''; end

