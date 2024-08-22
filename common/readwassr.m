function [B0map,hdr,pars,posimages,negimages,roimask,zspecint,ppmint,pathname] = readwassr(pathname,B0thresh,roimask,dofilt,swversion)

if nargin<2 || isempty(B0thresh), B0thresh = 1.0; end

if nargin<1 || isempty(pathname)
    [hdr,images,dicomhdr,pathname] = readdicomfiles2d;
    disp(pathname);
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname);
end
if nargin<4, dofilt = false; end

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
    if nargin<5 || isempty(swversion)
        hdr.swversion = 'VE11';
    else
        hdr.swversion = swversion;
    end
end

nx = size(images,1); ny = size(images,2);
nfreq = hdr.reps/2;
nz = size(images,3)/hdr.reps;
if nz>1
    images = reshape(images,[nx,ny,nz,nfreq*2]); 
end

if nargin<3 || isempty(roimask)
    if nz>1
        roimask = zeros(nx,ny,nz);
        for jj=1:nz
            roimask(:,:,jj) = images(:,:,jj,end)>mean(col(images(:,:,jj,end)))/2;
        end
    else
        roimask = images(:,:,end)>mean(col(images(:,:,end)))/2; % last image should have high signal
    end
end


seqname = hdr.SeqName;
swversion = hdr.swversion;
preptypeInd = []; cestz2value = []; wassrvalue = [];
if contains(swversion,'XA')
    if contains(seqname,'prepFLASH')
        preptypeInd = 2;
        wipppmindex = 2;
        cestz2value = 3;
        wassrvalue = 4;
        WASSRpw = hdr.WIPlong(13);
        WASSRpw1 = hdr.WIPlong(14);
        WASSRb1 = hdr.WIPlong(16);
        WASSRdc = hdr.WIPdbl(1);  
    else
        error('unknown XA sequence');
    end
elseif ( strfind( swversion, 'VD13A' ) )
    wipppmindex = 3;
%     prepmodeindex = 1;
    cestz2value = 5;
    WASSRpw = hdr.WIPlong(11);
    WASSRpw1 = hdr.WIPlong(13);
    WASSRb1 = hdr.WIPlong(14);
    WASSRdc = hdr.WIPdbl(1);
elseif ( contains( swversion, 'VD13D' ) || contains(swversion,'VE11') || contains(swversion,'VE12') )
    seqname = strsplit(seqname, '\');
    seqname = seqname{end};
    
    if strcmp(seqname,'prep_moco')
        wipppmindex = 2;
%         prepmodeindex = 2;
        cestz2value = 3;
        WASSRpw = 0.001*hdr.WIPlong(13);
        WASSRpw1 = 0.001*hdr.WIPlong(14);
        WASSRb1 = hdr.WIPlong(16);
        WASSRdc = hdr.WIPdbl(1);
    elseif strcmp(seqname,'prep_moco_iB0') || strcmp(seqname,'prep_moco_sB0')
        wipppmindex = 3; % iB0
%         prepmodeindex = 2;
        cestz2value = 3;
        WASSRpw = 0.001*hdr.WIPlong(15); % iB0
        WASSRpw1 = 0.001*hdr.WIPlong(16); % iB0
        WASSRb1 = hdr.WIPlong(18); % iB0
        WASSRdc = hdr.WIPdbl(1);
    elseif strcmp(seqname,'prep_tfl')
        wipppmindex = 2;
        prepmodeindex = 2;
        cestzvalue = 1e5;
        WASSRpw = hdr.WIPlong(13);
        WASSRpw1 = hdr.WIPlong(14);
        WASSRb1 = hdr.WIPlong(16);
        WASSRdc = hdr.WIPdbl(1);  
    elseif strcmp(seqname,'prepCV') || strcmp(seqname,'prep_tfl_nus')
        wipppmindex = 2;
        WASSRpw = hdr.WIPlong(14);
        WASSRpw1 = hdr.WIPlong(15);
        WASSRb1 = hdr.WIPlong(17);
        WASSRdc = hdr.WIPdbl(1);
	elseif strcmp(seqname,'prep_tfl_FatSat')
        wipppmindex = 2;
        WASSRpw = hdr.WIPlong(11);
        WASSRpw1 = hdr.WIPlong(12);
        WASSRb1 = hdr.WIPlong(14);
        WASSRdc = hdr.WIPdbl(1);										
    else
        wipppmindex = input(' Type in the value for wipppmindex : ');
%         prepmodeindex = input(' Type in the value for prepmodeindex : ');
        WASSRb1 = input(' Type in the value for CEST b1 in Hz : ');
        WASSRpw = input(' Type in the value for CEST Pulse duration in ms : ');
        WASSRpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
        WASSRdc = input(' Type in the value for CEST dutycycle in % : ');
    end
else
    if (strfind(seqname,'prep_moco'))
        wipppmindex = 2;
%         prepmodeindex = 2;
        cestz2value = 3;
        WASSRpw = 0.001*hdr.WIPlong(13);
        WASSRpw1 = 0.001*hdr.WIPlong(14);
        WASSRb1 = hdr.WIPlong(16);
        WASSRdc = hdr.WIPdbl(1);
    elseif (strfind(seqname,'prep_cv'))
        wipppmindex = 2;
%         prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            WASSRpw = hdr.WIPlong(10);
            WASSRb1 = hdr.WIPlong(13);
            WASSRpw1 = hdr.WIPlong(11);
            WASSRdc = hdr.WIPdbl(1);
        else
            WASSRb1 = input(' Type in the value for CEST b1 in Hz : ');
            WASSRpw = input(' Type in the value for CEST Pulse duration in ms : ');
            WASSRpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            WASSRdc = input(' Type in the value for CEST dutycycle in % : ');
        end
    elseif (strfind(seqname,'prep_flash'))
%         prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            WASSRpw = hdr.WIPlong(11);
            WASSRb1 = hdr.WIPlong(12);
            WASSRpw1 = hdr.WIPlong(15);
            WASSRdc = hdr.WIPdbl(1);
            wipppmindex = 3;
        else
            WASSRb1 = input(' Type in the value for CEST b1 in Hz : ');
            WASSRpw = input(' Type in the value for CEST Pulse duration in ms : ');
            WASSRpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            WASSRdc = input(' Type in the value for CEST dutycycle in % : ');
            wipppmindex = 5;
        end
    end
end

if ~isempty(preptypeInd) && hdr.WIPlong(preptypeInd)==wassrvalue && hdr.reps==9 % for fixed WASSR scans

    ppmstep  = hdr.WIPdbl(wipppmindex+3);  % step
    WASSRppmlist = 0:ppmstep:4*ppmstep;
    if nz>1
        posimages = images(:,:,:,end-4:end);
        negimages = flip(images(:,:,:,1:5),ndims(images));
        B0map = zeros(nx,ny,nz);
        for jj=1:nz
            [tempB0,tempzspec,ppmint] = mexwassr_b0map2d_extoutNW(WASSRppmlist,squeeze(posimages(:,:,jj,:)),squeeze(negimages(:,:,jj,:)),double(squeeze(roimask(:,:,jj))));
            tempB0 = anisodiff(tempB0,20,50,0.03,1);
            B0map(:,:,jj) = tempB0;
            zspecint(:,:,jj,:) = tempzspec; % FIX THIS
        end
    else
        posimages = images(:,:,end-4:end);
        negimages = flip(images(:,:,1:5),ndims(images));
        % B0map = mexwassr_b0map2d(WASSRppmlist, posimages, negimages, double(roimask)); % crashes without appropriate roimask
        [B0map, zspecint, ppmint] = mexwassr_b0map2d_extoutNW(WASSRppmlist, posimages, negimages, double(roimask)); % crashes without appropriate roimask
        B0map = anisodiff(B0map,20,50,0.03,1);
    end

else

    ppmbegin = hdr.WIPdbl(wipppmindex+1);  % begin
    ppmend   = hdr.WIPdbl(wipppmindex+2);  % end
    ppmstep  = hdr.WIPdbl(wipppmindex+3);  % step
    if (ppmbegin > ppmend)
        ppmstep = -ppmstep;
    end
    WASSRppmlist = ppmbegin:ppmstep:ppmend;
    if nz>1
        posimages = images(:,:,:,1:2:end);
        negimages = images(:,:,:,2:2:end);
        B0map = zeros(nx,ny,nz);
        for jj=1:nz
            [tempB0,tempzspec,ppmint] = mexwassr_b0map2d_extoutNW(WASSRppmlist,squeeze(posimages(:,:,jj,:)),squeeze(negimages(:,:,jj,:)),double(squeeze(roimask(:,:,jj))));
            tempB0 = anisodiff(tempB0,20,50,0.03,1);
            B0map(:,:,jj) = tempB0;
            zspecint(:,:,jj,:) = tempzspec; % FIX THIS
        end
    else
        posimages = images(:,:,1:2:end);
        negimages = images(:,:,2:2:end);
        % B0map = mexwassr_b0map2d(WASSRppmlist, posimages, negimages, double(roimask)); % crashes without appropriate roimask
        [B0map, zspecint, ppmint] = mexwassr_b0map2d_extoutNW(WASSRppmlist, posimages, negimages, double(roimask)); % crashes without appropriate roimask
        B0map = anisodiff(B0map,20,50,0.03,1);
    end

end

B0map( abs(B0map) > B0thresh ) = 0;
roimask( abs(B0map) > B0thresh ) = 0;

pars.ppmlist = WASSRppmlist;
pars.b1 = WASSRb1;
pars.dc = WASSRdc;
pars.pw = WASSRpw;
pars.pw1 = WASSRpw1;

