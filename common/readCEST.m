function [posimages, negimages, hdr, pars, pathname1] = readCEST(pathname1,dofilt,swversion)

if nargin<1 || isempty(pathname1)
    [hdr,images,dicomhdr,pathname1] = readdicomfiles2d;
elseif iscell(pathname1) % check this
    images = [];
%     pathname1 = NWsortscans(pathname1);
    for jj=1:length(pathname1)
        [hdr,temp,dicomhdr] = readdicomfiles2d(pathname1{jj});
        images = cat(3,images,temp);
    end
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname1);
end
if nargin<2, dofilt = false; end

if dofilt
    si = size(images);
    if length(dofilt)==3
        filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
    elseif length(dofilt)==6
        filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
    else
        filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
    end
    images = abs(fft2c(ifft2c(images).*filt));
end

if ~isfield(hdr,'swversion')
    if nargin<3 || isempty(swversion)
        hdr.swversion = 'VE11';
    else
        hdr.swversion = swversion;
    end
end
% 
% flag_mismatch = false;
% if ~isequal(hdr.reps,size(images,3))
%     nfreq = size(images,3)/2;
%     nfiles = size(images,3)/2;
%     warndlg('Reps in header do not match image size')
%     flag_mismatch = true;
% else
%     nfreq = hdr.reps/2;
%     nfiles = hdr.reps;
% end

is3d = false;
if ~isequal(hdr.reps,size(images,3))
    images = reshape(images,size(images,1),size(images,2),[],hdr.reps);
    is3d = true;
end
                    
seqname = hdr.SeqName;
swversion = hdr.swversion;
shotTRms = [];			  
if contains(swversion,'XA')
    if contains(seqname,'prepFLASH')
        wipppmindex = 2;
        CESTpw = hdr.WIPlong(13);
        CESTpw1 = hdr.WIPlong(14);
        CESTb1 = hdr.WIPlong(16);
        CESTdc = hdr.WIPdbl(1);  
        reqreps = hdr.WIPlong(9);
        mocoreps = hdr.WIPlong(11);
		shotTRms = hdr.WIPlong(8); % add this for other baselines													 
    else
        error('unknown XA sequence');
    end
elseif ( strfind( swversion, 'VD13A' ) )
    if (~strfind(seqname,'cest_flash'))
        warning('not cest_flash sequence')
        disp(seqname)
    end
    wipppmindex = 3;
    prepmodeindex = 1;
    cestz2value = 5;
    cestzvalue = 1e5; % arbitrarily large
    CESTpw = hdr.WIPlong(11);
    CESTpw1 = hdr.WIPlong(12);
    CESTb1 = hdr.WIPlong(14);
    CESTdc = hdr.WIPdbl(1);
    reqreps = hdr.WIPlong(8);
    mocoreps = 1;
elseif ( contains( swversion, 'VD13D' ) || contains(swversion,'VE11') || contains(swversion,'VE12') )
    seqname = strsplit(seqname, '\');
    seqname = seqname{end};
    if strcmp(seqname,'prep_moco')
        wipppmindex = 2;
        prepmodeindex = 2;
        cestz2value = 3;
        cestzvalue = 11; % NW
        CESTpw = 0.001*hdr.WIPlong(13);
        CESTpw1 = 0.001*hdr.WIPlong(14);
        CESTb1 = hdr.WIPlong(16);
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(9); % NW
        mocoreps = hdr.WIPlong(11); % NW
    elseif strcmp(seqname,'prep_moco_iB0') || strcmp(seqname,'prep_moco_sB0') 
        wipppmindex = 3; % iB0
        prepmodeindex = 2;
        cestz2value = 3;
        cestzvalue = 11; % NW
        CESTpw = 0.001*hdr.WIPlong(15); % iB0
        CESTpw1 = 0.001*hdr.WIPlong(16); % iB0
        CESTb1 = hdr.WIPlong(18); % iB0
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(10); % NW iB0
        mocoreps = hdr.WIPlong(12); % NW iB0
    elseif strcmp(seqname,'prep_tfl') || strcmp(seqname, 'prep_tfl_latest')
        wipppmindex = 2;
        prepmodeindex = 2;
        cestzvalue = 1e5;
        CESTpw = hdr.WIPlong(13);
        CESTpw1 = hdr.WIPlong(14);
        CESTb1 = hdr.WIPlong(16);
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(9);
        mocoreps = 1;  
        shotTRms = hdr.WIPlong(8); % add this for other baselines
    elseif strcmp(seqname,'prepCV') || strcmp(seqname,'prep_tfl_nus')
        wipppmindex = 2;
        prepmodeindex = 2;
        cestzvalue = 1e5;
        CESTpw = hdr.WIPlong(14);
        CESTpw1 = hdr.WIPlong(15);
        CESTb1 = hdr.WIPlong(17);
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(10);
        mocoreps = hdr.WIPlong(11);
    elseif strcmp(seqname,'prep_tfl_FatSat')
        wipppmindex = 2;
        prepmodeindex = nan;
        cestzvalue = 1e5;
        CESTpw = hdr.WIPlong(11);
        CESTpw1 = hdr.WIPlong(12);
        CESTb1 = hdr.WIPlong(14);
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(10);
        mocoreps = 1;  
    else
        wipppmindex = input(' Type in the value for wipppmindex : ');
        prepmodeindex = input(' Type in the value for prepmodeindex : ');
        cestzvalue = input(' Type in CESTz value : ');
        CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
        CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
        CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
        CESTdc = input(' Type in the value for CEST dutycycle in % : ');
    end
else
    cestzvalue = 1e5; % arbitrarily large
    if (strfind(seqname,'prep_moco'))
        wipppmindex = 2;
        prepmodeindex = 2;
        cestz2value = 3;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(10);
            CESTb1 = hdr.WIPlong(13);
            CESTpw1 = hdr.WIPlong(11);
            CESTdc = hdr.WIPdbl(1);
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
        end
    elseif (strfind(seqname,'prep_cv'))
        wipppmindex = 2;
        prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(10);
            CESTb1 = hdr.WIPlong(13);
            CESTpw1 = hdr.WIPlong(11);
            CESTdc = hdr.WIPdbl(1);
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
        end
    elseif (strfind(seqname,'prep_flash'))
        prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(11);
            CESTb1 = hdr.WIPlong(12);
            CESTpw1 = hdr.WIPlong(15);
            CESTdc = hdr.WIPdbl(1);
            wipppmindex = 3;
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
            wipppmindex = 5;
        end
    end
end
if ~exist('reqreps','var')
    warning('Required reps and moco reps were not read in')
    reqreps = 1;
    mocoreps = 1;
end
if is3d
    posimages = images(:,:,:,1:2:end);
    negimages = images(:,:,:,2:2:end);
else
    posimages = images(:,:,1:2:end);
    negimages = images(:,:,2:2:end);
end
ppmbegin = hdr.WIPdbl(wipppmindex+1);  % begin
ppmend   = hdr.WIPdbl(wipppmindex+2);  % end
ppmstep  = hdr.WIPdbl(wipppmindex+3);  % step
if (ppmbegin > ppmend)
    ppmstep = -ppmstep;
end
CESTppmlist = ppmbegin:ppmstep:ppmend;
% if hdr.WIPlong(prepmodeindex) == cestzvalue
%     DSimage = negimages(:,:,1);
%     posimages = posimages(:,:,2:2:end);
%     negimages = negimages(:,:,2:2:end);
% end
% if flag_mismatch
%     reqreps = size(posimages,3)/length(CESTppmlist);
%     mocoreps = 1;
% end

pars.CESTpw = CESTpw;
pars.CESTb1 = CESTb1;
pars.CESTpw1 = CESTpw1;
pars.CESTdc = CESTdc;
pars.reqreps = reqreps;
pars.mocoreps = mocoreps;
pars.CESTppmlist = CESTppmlist;
pars.shotTRms = shotTRms; 
pars.is3d = is3d;

end