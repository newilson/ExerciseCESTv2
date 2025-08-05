function [B1map,pars,hdr,alphaVec,images,B1fit,pathname1] = readB1reps(pathname1,dofilt,roimask)

nInterp = 50;

if nargin<1 || isempty(pathname1)
    [hdr,images,dicomhdr,pathname1] = readdicomfiles2d;
elseif isstruct(pathname1)
    hdr = pathname1.hdr;
    images = pathname1.images;
    if isfield(pathname1,'dicomhdr')
        dicomhdr = pathname1.dicomhdr;
    else
        dicomhdr = [];
    end
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname1);
end
if nargin<2, dofilt = false; end

if hdr.WIPlong(2)~=16
    warning('wrong sequence? double check')
end

reqreps = hdr.WIPlong(9);
actreps = hdr.WIPlong(10);

if actreps~=hdr.reps
    warning('repetition calculations do not match')
end

doLinFit = false;
% if actreps==reqreps+1
%     doLinFit = true;
% end

nx = size(images,1); ny = size(images,2);
nz = size(images,3)/hdr.reps;

if nargin<3 || isempty(roimask)
    % take MASK from first image which will have highest signal
    if nz>1
        roimask = ones(nx,ny,nz);
        for ii=1:nz
            roimask(:,:,ii) = images(:,:,ii,1) > mean(col(images(:,:,ii,1)))/2;
        end
    else
        roimask = images(:,:,1) > mean(col(images(:,:,1)))/2; 
    end
end


if dofilt
    si = size(images);
    if length(dofilt)==3
        if dofilt(1)~=si(1)
            dofilt(1) = si(1);
        end
        filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 si(3)]);
    elseif length(dofilt)==6
        if dofilt(1)~=si(1)
            dofilt(1) = si(1);
        end
        if dofilt(2)~=si(1)
            dofilt(2) = si(1);
        end
        filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 si(3)]);
    else
        filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
    end
    images = abs(fft2c(ifft2c(images).*filt));
end

% scale the images
if nz>1
    images_sc = images ./ (eps+repmat(images(:,:,:,1),[1 1 1 size(images,ndims(images))]));
else
    images_sc = images ./ (eps+repmat(images(:,:,1),[1 1 size(images,ndims(images))]));
end

alphaStart = hdr.WIPdbl(2)*pi/180;
alphaEnd = hdr.WIPdbl(7)*pi/180;
alphaVec = linspace(alphaStart,alphaEnd,reqreps);
if actreps==reqreps+1
    alphaVec = cat(2,0,alphaVec); % add first point
end
alphaVec = alphaVec(:);
alphaInterp = linspace(alphaVec(1),alphaVec(end), nInterp);

if doLinFit

    % Uses same model as nonlinear fit except that amplitude is known from the
    % signal when alpha = 0. Then acos(Sig/Amp) = B1rel x alpha where B1rel

    ind = find(alphaVec/1.5 > pi/2); % don't expect relative B1 > 1.5
    if isempty(ind)
        p = NWpolyfitim(1,alphaVec,acos(images_sc));
        if nz>1
            B1map = p(:,:,:,1);
        else
            B1map = p(:,:,1);
        end
        B1fit = NWpolyvalim(p,alphaInterp);
    else
        imagesFlip = images_sc;        
        for ii=1:length(ind)
            if nz>1
                imagesFlip(:,:,:,ind+ii-1:end) = -images_sc(:,:,:,ind+ii-1:end);
            else
                imagesFlip(:,:,ind+ii-1:end) = -images_sc(:,:,ind+ii-1:end);
            end
            pTemp = NWpolyfitim(1,alphaVec,acos(imagesFlip));


        end
    end

else

    % pars:
    %    ampl
    %    relB1
    B1func = @(par,alpha) par(1) * cos(par(2) * alpha);
    lb = [0.7 0];
    ub = [2 3];
    par0 = [1.3 1];
    options = optimset('Algorithm','trust-region-reflective', 'MaxIter', 100, 'MaxFunEvals', 200, 'TolFun',1e-3,'TolX',1e-6,'display', 'off');


    B1fit = zeros([size(roimask) nInterp]);
    B1map = zeros(size(roimask));
    pars = zeros([size(roimask) length(par0)]);

    for ii=1:nx*ny*nz
        if nz>1
            [x,y,z] = ind2sub([nx,ny,nz],ii);
        else
            [x,y] = ind2sub([nx,ny],ii);
            z = 1;
        end

        if roimask(x,y,z)==1
            %do fit
            if nz>1
                tempVec = squeeze(images_sc(x,y,z,:));
            else
                tempVec = squeeze(images_sc(x,y,:));
            end
            tempVec = tempVec(:);
            [~,minInd] = min(tempVec);

            % Min Signal is negative
            tempVec1 = tempVec; tempVec1(minInd:end) = -tempVec(minInd:end);
            par1 = lsqcurvefit(B1func,par0,alphaVec,tempVec1,lb,ub,options); % uses previous iteration fitted parameters as next voxel initial guess
            fitVec1 = B1func(par1,alphaVec);
            gof1 = goodnessOfFit_(fitVec1,tempVec1,'MSE');

            % Min Signal is positive
            tempVec2 = -tempVec; tempVec2(1:minInd) = tempVec(1:minInd);
            par2 = lsqcurvefit(B1func,par0,alphaVec,tempVec2,lb,ub,options);
            fitVec2 = B1func(par2,alphaVec);
            gof2 = goodnessOfFit_(fitVec2,tempVec2,'MSE');

            % Nothing is negative - can be different from 2nd case for
            % noisy data
            tempVec3 = tempVec;
            par3 = lsqcurvefit(B1func,par0,alphaVec,tempVec3,lb,ub,options);
            fitVec3 = B1func(par3,alphaVec);
            gof3 = goodnessOfFit_(fitVec3,tempVec3,'MSE');

            minGof = min([gof1 gof2 gof3]);
            if minGof==gof1
                par = par1;
            elseif minGof==gof2
                par = par2;
            else
                par = par3;
            end

            fitVec = B1func(par,alphaInterp(:));
            if nz>1
                B1fit(x,y,z,:) = fitVec * images(x,y,z,1);
                B1map(x,y,z) = par(2);
                pars(x,y,z,:) = par;
            else
                B1fit(x,y,:) = fitVec * images(x,y,1);
                B1map(x,y) = par(2);
                pars(x,y,:) = par;
            end
        end
    end

end

if nz>1
    for jj=1:nz
        B1map(:,:,jj) = anisodiff(B1map(:,:,jj),20,50,0.03,1);
    end
else
    B1map = anisodiff(B1map,20,50,0.03,1);
end

