function [T2star,pars,hdr,TEms,images,fitImages,gof,pathname] = readT2star(pathname,dofilt)

if nargin<2 || isempty(dofilt), dofilt = false; end

if nargin<1 || isempty(pathname)
    [hdr,images,dicomhdr,pathname] = readdicomfiles2d;
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname);
end

nz = size(images,3)/hdr.reps/hdr.contrasts;

TEms = hdr.TEms(1:hdr.contrasts);

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

% reshape
nx = size(images,1); ny = size(images,2);
images = reshape(images,nx,ny,nz,hdr.contrasts,hdr.reps); % double check this if slices come before contrasts
images = permute(images,[1 2 3 5 4]);

% linear fit
logI = log(images);
poly_deg = 1;
pars = NWpolyfitim(poly_deg,TEms,logI);

% evaluation
logIfit = NWpolyvalim(pars,TEms);

T2star = -1./pars(:,:,:,:,1);
fitImages = exp(logIfit);
gof = goodnessOfFit_(fitImages(:),images(:),'NRMSE'); % note that gof is reported on the images, not log(images) which are actually fit
