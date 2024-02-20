% Function standalone_sort
% Author: Matthew A. Sochor - matthew.sochor@gmail.com
% Date: 3/28/2008
%
%  This function sorts a folder of dicoms (as you might
%  recieve from the COMQUEST dicom server or similar) into
%  the following folder structure:
%
%  Patient ID\Study Date\ImageSeries\ (file name)
%
%  Syntax:
%     standalone_sort( [path1] );
%     standalone_sort( [path1],[writepath1] );
%     standalone_sort( [path1],[writepath1],[copyflag] );
%     standalone_sort( [path1],[writepath1],[copyflag],[v2check] );
%
%  Input:
%     path1: (optional) folder containing dicom images
%     writepath1: (optional) folder to write new folder structure into
%     copyflag: (optional) 1: Keep original file and copy to new location
%                          0: Move file to new location
%     v2check: (optional)  1: Check if file is a '.v2' ending DICOM.
%                          0: Copy any file in the folder.  NOTE: this can
%                          cause errors if the folder has non-DICOM images
%
%   Thoughts:  To improve robustness, there should be a dcmcheck, I just
%   didn't do it.  Also maybe make it a switch/case style to be less
%   rigorous on the input order.

% NW 2/6/2017


function outdir = standalone_sort(varargin)
if nargin < 1
    path1 = uigetdir('Open directory containing DICOM images');
    writepath1 = path1;
    copyflag = 1;
    v2check = 0;
elseif nargin == 1
    path1 = varargin{1};
    writepath1 = path1;
    copyflag = 1;
    v2check = 0;
elseif nargin == 2
    path1 = varargin{1};
    writepath1 = varargin{2};
    copyflag = 1;
    v2check = 0;
elseif nargin == 3
    path1 = varargin{1};
    writepath1 = varargin{2};
    copyflag = varargin{3};
    v2check = 0;
elseif nargin == 4
    path1 = varargin{1};
    writepath1 = varargin{2};
    copyflag = varargin{3};
    v2check = varargin{4};
elseif nargin > 4
    disp('Only four or less input arguments allowed.  Returning...');
    return;
end

tic
ext     = '.dcm';		% the output filename extension
files = dir([path1 filesep '**' filesep '*']);
namestr = [];
for j=1:length(files)
    if ~files(j).isdir
%     fullfilen = [path1 filesep files(j).name];
    fullfilen = fullfile(files(j).folder,files(j).name);
    if ( ~isdir(fullfilen)  && isdicom (fullfilen) )
        info = dicominfo(fullfilen);
        if (isempty(namestr))
            if( isfield(info.PatientName,'GivenName') && isfield(info.PatientName,'FamilyName') )
                namestr = [info.PatientName.GivenName info.PatientName.FamilyName];
            elseif ( isfield(info.PatientName,'GivenName') )
                namestr = info.PatientName.GivenName;
            elseif ( isfield(info.PatientName,'FamilyName') )
                namestr = info.PatientName.FamilyName;
            else
                namestr = 'ANONYMOUS';
            end
            namestr = regexprep(namestr,' ','');
            outdirtop = [writepath1 filesep namestr '-' info.StudyDate];
            if (~isdir(outdirtop))
                [success ,message, messageid] = mkdir(outdirtop);
            end
        end
        exnum   = str2num(info.StudyID);
        if (isempty(exnum))
            exnum = 1;
        end
        outdir = fullfile( outdirtop,sprintf('E%05i',exnum*100) );
        if ~isfield(info,'SeriesDescription') % anonymized? or xa? not sure (NW)
            pvthdr = char(info.SharedFunctionalGroupsSequence.Item_1.Private_0021_10fe.Item_1.Private_0021_1019)'; % XA enhanced dicom
            ind = strfind(pvthdr,'tProtocolName');
            ind2 = strfind(pvthdr(ind(1):end),'""');
            outstr = pvthdr(1+ind+ind2(1):ind+ind2(2)-2);
        else
            outstr = info.SeriesDescription;
        end
        outstr = regexprep(outstr,'<','');
        outstr = regexprep(outstr,'>','');
%         outdir = fullfile( outdir, sprintf('S%05i-%s',info.SeriesNumber,outstr) );
        outdir = fullfile(outdir, sprintf('%s_%04i',outstr,info.SeriesNumber)); % NW
        if (~isdir(outdir))
            [success ,message, messageid] = mkdir(outdir);
        end
        
%         if (isfield(info,'EchoNumber'))
%             outname = sprintf('E%05iS%05iX%1iI%05i',exnum*100,info.SeriesNumber,info.EchoNumber,info.InstanceNumber);
%         else
%             outname = sprintf('E%05iS%05iX1I%05i',exnum*100,info.SeriesNumber,info.InstanceNumber);
%         end
        outname = sprintf('E%03iS%03iX1I%03i',exnum,info.SeriesNumber,info.InstanceNumber); % NW

        
        if copyflag
            [status,message,messageid] = copyfile(fullfile(files(j).folder,files(j).name),[outdir filesep outname ext]);
        else
            [status,message,messageid] = movefile(fullfile(files(j).folder,files(j).name),[outdir filesep outname ext]);
        end
    end
    end
end

toc
