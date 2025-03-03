function [refimage,hdr,dicomhdr,refimageb,mask,biasscale,pathname] = readref(pathname,dofilt)

if nargin<1 || isempty(pathname)
    [hdr, refimage,dicomhdr,pathname] = readdicomfiles2d;
elseif isstruct(pathname)
    hdr = pathname.hdr;
    refimage = pathname.refimage;
    if isfield(pathname,'dicomhdr')
        dicomhdr = pathname.dicomhdr;
    else
        dicomhdr = [];
    end
else
    [hdr,refimage,dicomhdr] = readdicomfiles2d(pathname);
end
if nargin<2, dofilt = false; end

if dofilt
    si = size(refimage);
    if length(dofilt)==3
        filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
    elseif length(dofilt)==6
        filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
    else
        filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
    end
    refimage = abs(fft2c(ifft2c(refimage).*filt));
end

[H,W,Z] = size(refimage);

if 1 % always filter
    for ii=1:Z
        refimage(:,:,ii) = anisodiff(refimage(:,:,ii),10,50,0.03,1);
    end
end
refimage_out = refimage;
refimage = 2000.0 * refimage/max(refimage(:));
noise = refimage(refimage < 100); % was 50
noisestd = std(noise(:));
threshold = 4 * noisestd; % was 8
biasimages2 = 0*refimage;
for ii=1:Z
    biasimages2(:,:,ii) = anisodiff(refimage(:,:,ii),10,10000,0.25,1);
end
% mask = zeros(H,W);
mask = ones(size(refimage));
im1 = refimage;
maxval1 = max(im1(:));
mask(im1 < threshold) = 0.0;
im = biasimages2;
maxval = max(im(:));
im = maxval ./ (im+0.01);
im (mask == 0 ) = 1.0;
im (im < 0) = 1.0;
im (im > 30) = 1.0; % was 50
im1 =  (refimage .* mask) .* im;
im(im1 > maxval1) = im(im1 > maxval1) .* (maxval1 ./ (im1(im1 > maxval1)));  % Division to obtain biasscale
biasscale = im;
refimageb = refimage .* biasscale;
refimageb = 2000.0 * refimageb/max(refimageb(:));
refimage = refimage_out;
