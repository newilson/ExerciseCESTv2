function [CESTmap,B0corrp,B0corrn,inter] = getCESTFromScans(pars,CESTp,CESTn,B0map,B1map,roimask,opt,poly_deg)

if nargin<6 || isempty(roimask), roimask = ones(size(B1map)); end
if nargin<7, opt = 2; end
if nargin<8, poly_deg = 2; end

ppm = pars.CESTppm;
ppmlist = pars.CESTppmlist;
reps = size(CESTp,3)/length(ppmlist);

CESTmap = zeros(size(CESTp,1),size(CESTp,2),reps);
for ii=1:reps
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
            im1 = cat(3,flip(nimg,3),pimg);
            ax1 = cat(2,-flip(ppmlist),ppmlist);
            im2 = cat(3,flip(nimginterp,3),pimginterp);
            ax2 = cat(2,-flip(ppminterp),ppminterp);
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
    B1corrmap = roimask .* (B0corrmap ./ (B1map+eps)); % linear B1 correction
    B1corrmap(~isfinite(B1corrmap)) = 0;
%     B1corrmap = anisodiff(B1corrmap,20,50,0.03,1);
    
    CESTmap(:,:,ii) = B1corrmap;
end 


