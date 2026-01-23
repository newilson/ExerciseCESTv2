function [B1map,pars,hdr,alphaVec,images,B1fit,alphaInterp,pathname1] = readB1reps(pathname1,dofilt,roimask,doLinFit)

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
if nargin<2 || isempty(dofilt), dofilt = false; end

if hdr.WIPlong(2)~=16
    warning('wrong sequence? double check')
end

reqreps = hdr.WIPlong(9);
actreps = hdr.WIPlong(10);

if actreps~=hdr.reps
    warning('repetition calculations do not match')
end

% if actreps==reqreps+1
%     doLinFit = true;
% end

nx = size(images,1); ny = size(images,2);
nz = size(images,3)/hdr.reps;

if nargin<4 || isempty(doLinFit)
    if nz>1
        doLinFit = true;
    else
        doLinFit = false;
    end
end

if nz>1
    images = reshape(images,nx,ny,nz,hdr.reps);
end

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

    images_sc = min(images_sc,1); % nothing bigger than 1 - could be a problem in noisy voxels, doesn't do anything for high SNR voxels -> forces acos to be real

    ind = find(alphaVec*1.5 > pi/2, 1); % possible negative signals - don't expect relative B1 > 1.5
    if isempty(ind) % expects no negative signal intensities
        p = NWpolyfitim(1, alphaVec, acos(images_sc));
        if nz>1
            B1map = p(:,:,:,1);
        else
            B1map = p(:,:,1);
        end
        B1fit = cos( NWpolyvalim(p, alphaInterp) );
    else
        % Sign ambiguity: try different crossing points, pick best per voxel
        % Baseline: no sign correction
        [p_best, S_best] = NWpolyfitim(1, alphaVec, acos(images_sc));

        % Try negating from each possible crossing point onward
        for jj = ind:length(alphaVec)
            imagesFlip = images_sc;
            if nz>1
                imagesFlip(:,:,:,jj:end) = -images_sc(:,:,:,jj:end);
            else
                imagesFlip(:,:,jj:end) = -images_sc(:,:,jj:end);
            end
            [pTemp, STemp] = NWpolyfitim(1, alphaVec, acos(max(imagesFlip, -1)));

            % Update voxels where this fit has lower residual
            better = STemp.normr < S_best.normr;
            if nz>1
                for kk = 1:size(p_best, 4)
                    temp = p_best(:,:,:,kk);
                    tempNew = pTemp(:,:,:,kk);
                    temp(better) = tempNew(better);
                    p_best(:,:,:,kk) = temp;
                end
            else
                for kk = 1:size(p_best, 3)
                    temp = p_best(:,:,kk);
                    tempNew = pTemp(:,:,kk);
                    temp(better) = tempNew(better);
                    p_best(:,:,kk) = temp;
                end
            end
            S_best.normr(better) = STemp.normr(better);
        end

        p = p_best;
        if nz>1
            B1map = p(:,:,:,1);
        else
            B1map = p(:,:,1);
        end
        B1fit = cos( NWpolyvalim(p, alphaInterp) );
    end

    % Scale B1fit by first image to get absolute signal
    if nz>1
        B1fit = B1fit .* repmat(images(:,:,:,1), [1 1 1 nInterp]);
    else
        B1fit = B1fit .* repmat(images(:,:,1), [1 1 nInterp]);
    end
    pars = p;

else

    % pars:
    %    ampl
    %    relB1
    B1func = @(par,alpha) par(1) * cos(par(2) * alpha);
    lb = [0.7 0];
    ub = [2 3];
    par0 = [1.3 1];
    npar = length(par0);
    options = optimset('Algorithm','trust-region-reflective', 'MaxIter', 100, 'MaxFunEvals', 200, 'TolFun',1e-3,'TolX',1e-6,'display', 'off');

    % Reshape for parallel-friendly linear indexing
    nVox = nx*ny*nz;
    images_sc_vec = reshape(images_sc, nVox, []);
    if nz>1
        images_first = reshape(images(:,:,:,1), nVox, 1);
    else
        images_first = reshape(images(:,:,1), nVox, 1);
    end
    roimask_vec = roimask(:);

    B1map_vec = zeros(nVox, 1);
    pars_vec = zeros(nVox, npar);
    B1fit_vec = zeros(nVox, nInterp);

    % Check for parallel computing toolbox
    doParFor = false;
    try
        checkParTool = canUseParallelPool;
    catch
        checkParTool = false;
    end
    if checkParTool
        if ~isempty(gcp('nocreate')) || nz>1
            doParFor = true;
        end
    end

    if doParFor
        parfor ii=1:nVox
            if roimask_vec(ii)==1
                tempVec = images_sc_vec(ii,:).';
                [~,minInd] = min(tempVec);

                % Min Signal is negative
                tempVec1 = tempVec; tempVec1(minInd:end) = -tempVec(minInd:end);
                par1 = lsqcurvefit(B1func,par0,alphaVec,tempVec1,lb,ub,options);
                fitVec1 = B1func(par1,alphaVec);
                gof1 = goodnessOfFit_(fitVec1,tempVec1,'MSE');

                % Min Signal is positive
                tempVec2 = -tempVec; tempVec2(1:minInd) = tempVec(1:minInd);
                par2 = lsqcurvefit(B1func,par0,alphaVec,tempVec2,lb,ub,options);
                fitVec2 = B1func(par2,alphaVec);
                gof2 = goodnessOfFit_(fitVec2,tempVec2,'MSE');

                % Nothing is negative
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

                B1fit_vec(ii,:) = B1func(par,alphaInterp(:)).' * images_first(ii);
                B1map_vec(ii) = par(2);
                pars_vec(ii,:) = par;
            end
        end
    else
        for ii=1:nVox
            if roimask_vec(ii)==1
                tempVec = images_sc_vec(ii,:).';
                [~,minInd] = min(tempVec);

                % Min Signal is negative
                tempVec1 = tempVec; tempVec1(minInd:end) = -tempVec(minInd:end);
                par1 = lsqcurvefit(B1func,par0,alphaVec,tempVec1,lb,ub,options);
                fitVec1 = B1func(par1,alphaVec);
                gof1 = goodnessOfFit_(fitVec1,tempVec1,'MSE');

                % Min Signal is positive
                tempVec2 = -tempVec; tempVec2(1:minInd) = tempVec(1:minInd);
                par2 = lsqcurvefit(B1func,par0,alphaVec,tempVec2,lb,ub,options);
                fitVec2 = B1func(par2,alphaVec);
                gof2 = goodnessOfFit_(fitVec2,tempVec2,'MSE');

                % Nothing is negative
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

                B1fit_vec(ii,:) = B1func(par,alphaInterp(:)).' * images_first(ii);
                B1map_vec(ii) = par(2);
                pars_vec(ii,:) = par;
            end
        end
    end

    % Reshape outputs back to image dimensions
    B1map = reshape(B1map_vec, size(roimask));
    pars = reshape(pars_vec, [size(roimask) npar]);
    B1fit = reshape(B1fit_vec, [size(roimask) nInterp]);

end

if nz>1
    for jj=1:nz
        B1map(:,:,jj) = anisodiff(B1map(:,:,jj),20,50,0.03,1);
    end
else
    B1map = anisodiff(B1map,20,50,0.03,1);
end

