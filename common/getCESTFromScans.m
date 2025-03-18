function [imgs,inter] = getCESTFromScans(pars,imgs,opt,poly_deg)
% B0map must be [x,y,z] or [x,y,z,off,1:2] where the 5th dimension index 1
% is the positive offsets, and index 2 is the negative offsets

if nargin<3, opt = 2; end
if nargin<4, poly_deg = 2; end
if ~isfield(pars,'b1corrMode') || isempty(pars.b1corrMode)
    pars.b1corrMode = 'linear';
end
if isfield(imgs,'CESTp'), CESTp = imgs.CESTp; else warning('cannot get CEST without CESTp images'); return; end
if isfield(imgs,'CESTn'), CESTn = imgs.CESTn; else warning('cannot get CEST without CESTn images'); return; end
si = size(CESTp);
if isfield(imgs,'B0map'), B0map = imgs.B0map; else warning('no B0map used'); B0map = zeros(si(1:end-1)); end
if isfield(imgs,'B1map'), B1map = imgs.B1map; else warning('no B1map used'); B1map = ones(si(1:end-1)); end
if isfield(imgs,'roimask'), roimask = imgs.roimask; else roimask = ones(si(1:end-1)); end
if isfield(imgs,'T1map'), T1map = imgs.T1map; else T1map = []; end



ppm = pars.CESTppm;
ppmlist = pars.CESTppmlist;
reps = size(CESTp,ndims(CESTp))/length(ppmlist);
if isfield(pars,'CESTb1') && ~isempty(pars.CESTb1)
    CESTb1 = pars.CESTb1;
else
    CESTb1 = 200; % Hz?
end


if ndims(CESTp)>3
    is3d = true;
    CESTmap = zeros(size(CESTp,1),size(CESTp,2),size(CESTp,3),reps);
else
    is3d = false;
    CESTmap = zeros(size(CESTp,1),size(CESTp,2),reps);
end
B0CorrP = 0*CESTmap; B0CorrN = 0*CESTmap;
inter = [];

% multiple repetitions
for ii=1:reps

    % B0 correction
    if is3d
        pimg = CESTp(:,:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
        nimg = CESTn(:,:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));

        if size(B0map,4)==1
            if opt>2
                warning('only OPTION 2 is supported, using that one')
                opt = 2;
            end
            [pout,nout,pfit,nfit] = B0correctedCESTimg_ASNW(abs(ppm),ppmlist,pimg,nimg,B0map,roimask,poly_deg,opt);
            ppminterp = linspace(min(ppmlist(:)),max(ppmlist(:)),30);
            pimginterp = NWpolyvalim(pfit,ppminterp);
            nimginterp = NWpolyvalim(nfit,ppminterp);
            inter.ppm = ppminterp;
            inter.pimg = pimginterp;
            inter.nimg = nimginterp;
            % im1 = cat(4,flip(nimg,4),pimg);
            % ax1 = cat(2,-flip(ppmlist),ppmlist);
            % im2 = cat(4,flip(nimginterp,4),pimginterp);
            % ax2 = cat(2,-flip(ppminterp),ppminterp);
            %             NWinteractiveCEST3d(im1,ax1,im2,ax2);
        else
            pB0 = B0map(:,:,:,:,1);
            nB0 = B0map(:,:,:,:,2);
            pB0 = pB0(:,:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
            nB0 = nB0(:,:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
            pout = zeros(size(CESTp,1),size(CESTp,2),size(CESTp,3));
            nout = zeros(size(CEStp,1),size(CESTp,2),size(CESTp,3));
            for ss=1:size(CESTp,3)
                [pout(:,:,ss), nout(:,:,ss)] = mexB0correctedCEST2dmoco(abs(ppm),ppmlist,squeeze(pimg(:,:,ss,:)),squeeze(nimg(:,:,ss,:)),squeeze(pB0(:,:,ss,:)),squeeze(nB0(:,:,ss,:)),double(roimask(:,:,ss)));
            end
        end

        for ss = 1:size(CESTp,3)
            B0corrp(:,:,ss) = anisodiff(pout(:,:,ss),20,50,0.03,1);
            B0corrn(:,:,ss) = anisodiff(nout(:,:,ss),20,50,0.03,1);
        end

    else
        pimg = CESTp(:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
        nimg = CESTn(:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));

        if size(B0map,3)==1
            if opt>2
                [pout,nout] = mexB0correctedCESTimglfit(abs(ppm),ppmlist,pimg,nimg,B0map,double(roimask));
            else
                [pout,nout,pfit,nfit] = B0correctedCEST2Dimg_ASNW(abs(ppm),ppmlist,pimg,nimg,B0map,roimask,poly_deg,opt);
                ppminterp = linspace(min(ppmlist(:)),max(ppmlist(:)),30);
                pimginterp = NWpolyvalim(pfit,ppminterp);
                nimginterp = NWpolyvalim(nfit,ppminterp);
                inter.ppm = ppminterp;
                inter.pimg = pimginterp;
                inter.nimg = nimginterp;
                % im1 = cat(3,flip(nimg,3),pimg);
                % ax1 = cat(2,-flip(ppmlist),ppmlist);
                % im2 = cat(3,flip(nimginterp,3),pimginterp);
                % ax2 = cat(2,-flip(ppminterp),ppminterp);
                %             NWinteractiveCEST(im1,ax1,im2,ax2);
            end
        else
            pB0 = B0map(:,:,:,1);
            nB0 = B0map(:,:,:,2);
            pB0 = pB0(:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
            nB0 = nB0(:,:,1+(ii-1)*length(ppmlist):ii*length(ppmlist));
            [pout,nout] = mexB0correctedCEST2dmoco(abs(ppm),ppmlist,pimg,nimg,pB0,nB0,double(roimask));
        end

        B0corrp = anisodiff(pout,20,50,0.03,1);
        B0corrn = anisodiff(nout,20,50,0.03,1);

    end

    % B1 correction and MTR calculation
    if strcmpi(pars.b1corrMode,'gluCEST_NMRB') && isempty(T1map)
        warning('no T1map - doing LINEAR B1 correction');
        pars.b1corrMode = 'linear';
    end

    if strcmpi(pars.b1corrMode,'gluCEST_NMRB') % nonlinear - from NMRB paper - B1 correct images, then take asymmetry

        if ppm-3>0.1
            warning('gluCEST corrections should be applied at 3 ppm')
        end

        % load from mat file
        in = load('gluCESTB1cal_NMRB.mat');
        NumMasks = in.NumMasks;
        Intervals = in.Intervals;
        B1FitPars = in.B1FitPars;

        % Pseudo code - don't loop over pixels, loop over intervals
        % binning
        T1index = 0*T1map;
        for aa=1:NumMasks
            T1index = T1index + (T1map>Intervals(1,aa));
        end
        T1index(T1index==0) = 1; % bin to first - maybe this should be masked out instead

        %b1 correct
        if is3d
            H = si(1); W = si(2); D = si(3);
        else
            H = si(1); W = si(2); D = 1;
        end
        for k = 1:D
            roimaskSl = roimask(:,:,k);
            for j = 1:W
                for i = 1:H
                    if (roimaskSl(i,j) > 0 && B1map(i,j,k)>eps ) % NW
                        relB1 = B1map(i,j,k);
                        maskIndex = T1index(i,j,k); %identifies T1value interval of that pixel
                        if (maskIndex)
                            parameters = B1FitPars(maskIndex,:);
                            Bn = parameters(1); Cn = parameters(2); Dn = parameters(3); En = parameters(4);
                            Bp = parameters(5); Cp = parameters(6); Dp = parameters(7); Ep = parameters(8);
                            x = CESTb1;
                            refN = En*(1+(Bn*x.^2 ./ (Cn*x.^2 +1)) - Dn*x.^2);
                            refP = Ep*(1+(Bp*x.^2 ./ (Cp*x.^2 +1)) - Dp*x.^2);
                            x = CESTb1* relB1;
                            valN = En*(1+(Bn*x.^2 ./ (Cn*x.^2 +1)) - Dn*x.^2);
                            valP = Ep*(1+(Bp*x.^2 ./ (Cp*x.^2 +1)) - Dp*x.^2);
                            B0corrp(i,j,k) = B0corrp(i,j,k) * refP/valP;
                            B0corrn(i,j,k) = B0corrn(i,j,k) * refN/valN;
                        else
                            B0corrp(i,j,k) = 0.0;
                            B0corrn(i,j,k) = 0.0;
                        end
                    else
                        B0corrp(i,j,k) = 0.0;
                        B0corrn(i,j,k) = 0.0;
                    end
                end
            end
        end

        %mtr asym
        if ppm>0
            if isfield(pars,'norm') && ~isempty(pars.norm) && strcmpi(pars.norm,'m0') && isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (B0corrn - B0corrp) ./ (pars.normimg + eps); % user supplied normalization image
            elseif isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (pars.normimg - B0corrp) ./ (pars.normimg + eps); % user supplied normalization image
            else
                B0corrmap = 100 * (B0corrn - B0corrp) ./ (B0corrn+eps); % default - negative normalization
                pars.norm = 'neg';
                pars.normimg = B0corrn;
            end
        else % interested in negative side
            if isfield(pars,'norm') && ~isempty(pars.norm) && strcmpi(pars.norm,'m0') && isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (B0corrp - B0corrn) ./ (pars.normimg + eps); % user supplied normalization image
            elseif isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (pars.normimg - B0corrn) ./ (pars.normimg + eps); % user supplied normalization image
            else
                B0corrmap = 100 * (B0corrp - B0corrn) ./ (B0corrp+eps); % default - negative normalization
                pars.norm = 'neg';
                pars.normimg = B0corrp;
            end
        end

        if is3d
            CESTmap(:,:,:,ii) = B0corrmap;
            B0CorrP(:,:,:,ii) = B0corrp;
            B0CorrN(:,:,:,ii) = B0corrn;
        else
            CESTmap(:,:,ii) = B0corrmap;
            B0CorrP(:,:,ii) = B0corrp;
            B0CorrN(:,:,ii) = B0corrn;
        end

    else % linear/none - take asymmetry first, then B1 correct contrast

        if ppm>0
            if isfield(pars,'norm') && ~isempty(pars.norm) && strcmpi(pars.norm,'m0') && isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (B0corrn - B0corrp) ./ (pars.normimg + eps); % user supplied normalization image
            elseif isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (pars.normimg - B0corrp) ./ (pars.normimg + eps); % user supplied normalization image
            else
                B0corrmap = 100 * (B0corrn - B0corrp) ./ (B0corrn+eps); % default - negative normalization
                pars.norm = 'neg';
                pars.normimg = B0corrn;
            end
        else % interested in negative side
            if isfield(pars,'norm') && ~isempty(pars.norm) && strcmpi(pars.norm,'m0') && isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (B0corrp - B0corrn) ./ (pars.normimg + eps); % user supplied normalization image
            elseif isfield(pars,'normimg') && isequal(size(pars.normimg),size(B0corrn))
                B0corrmap = 100 * (pars.normimg - B0corrn) ./ (pars.normimg + eps); % user supplied normalization image
            else
                B0corrmap = 100 * (B0corrp - B0corrn) ./ (B0corrp+eps); % default - negative normalization
                pars.norm = 'neg';
                pars.normimg = B0corrp;
            end
        end

        B0corrmap(~isfinite(B0corrmap)) = 0;
        if strcmpi(pars.b1corrMode,'none')
            B1corrmap = roimask .* B0corrmap;
        else
            B1corrmap = roimask .* (B0corrmap ./ (B1map+eps)); % linear B1 correction
        end
        B1corrmap(~isfinite(B1corrmap)) = 0;
        %     B1corrmap = anisodiff(B1corrmap,20,50,0.03,1);

        if is3d
            CESTmap(:,:,:,ii) = B1corrmap;
            B0CorrP(:,:,:,ii) = B0corrp;
            B0CorrN(:,:,:,ii) = B0corrn;
        else
            CESTmap(:,:,ii) = B1corrmap;
            B0CorrP(:,:,ii) = B0corrp;
            B0CorrN(:,:,ii) = B0corrn;
        end

    end

end

imgs.CESTmap = CESTmap;
imgs.B0CorrP = B0CorrP;
imgs.B0CorrN = B0CorrN;

