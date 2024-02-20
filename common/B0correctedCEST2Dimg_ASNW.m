%% B0correctedCEST2Dimg_AS.m functing created by @Anup, July 2017 
% New function for correcting CEST weighted images for B0 inhomogenity
% 'ppm' at which outpot image is required
%  ppmlist-- list of ppm values at which input posImgs  were acquired, 
%  posImg, images at ppmlist
%  NegImg, images at -ppmlist
%  B0 map,  mask
%  poly_deg  --defalt = 2, (it can be increased to 3 or 4 depending upon
%  nature of posimg intensit vs ppmlist curve shape.

function [pimgout,nimgout,pfit,nfit] = B0correctedCEST2Dimg_ASNW(ppmOI,ppmlist,posimages,negimages,B0map,mask,poly_deg,opt)

if nargin<8 || isempty(opt), opt = 2; end
if ndims(B0map)>2, opt = 2; end % only option for general sizes

[W, H, nuImg] =size(posimages);
pimgout = zeros(W, H);
nimgout = zeros(W, H);
% NW
pfit = zeros(W,H,poly_deg+1);
nfit = zeros(W,H,poly_deg+1);

ppmlist = squeeze(ppmlist)';

% NW added faster option
switch opt
    case 0 % original - offset frequencies corrected before polyfit
        for i = 1:W
            for j= 1:H
                if(mask(i, j)>0)
                    B0offset = B0map(i,j);
                    
                    %%%% If both pos and neg images are from low to high
                    ppmB0corrP = ppmlist-B0offset;
                    ppmB0corrN = ppmlist+B0offset;
                    %%%%%
                    
                    posarr = squeeze(posimages(i,j, :));
                    negarr = squeeze(negimages(i,j, :));
                    
                    ppfit =  polyfit(ppmB0corrP, posarr, poly_deg);
                    npfit =  polyfit(ppmB0corrN, negarr, poly_deg);
                    
                    % NW
                    pfit(i,j,:) = ppfit;
                    nfit(i,j,:) = npfit;
                    
                    posB0corr_ppm = polyval(ppfit, ppmOI);
                    negB0corr_ppm = polyval(npfit, ppmOI);
                    
                    pimgout(i,j) = posB0corr_ppm;
                    nimgout(i,j) = negB0corr_ppm;
                                        
                end
            end
        end
        
    case 1 % NW new - offset frequencies corrected in polyval (amounts to shifting the desired frequency)
        for i = 1:W
            for j= 1:H
                if(mask(i, j)>0)
                    B0offset = B0map(i,j);
                    
                    %%%% If both pos and neg images are from low to high
                    ppmB0corrP = ppmlist;
                    ppmB0corrN = ppmlist;
                    %%%%%
                    
                    posarr = squeeze(posimages(i,j, :));
                    negarr = squeeze(negimages(i,j, :));
                    
                    ppfit =  polyfit(ppmB0corrP, posarr, poly_deg);
                    npfit =  polyfit(ppmB0corrN, negarr, poly_deg);
                    
                    % NW
                    pfit(i,j,:) = ppfit;
                    nfit(i,j,:) = npfit;
                    
                    posB0corr_ppm = polyval(ppfit, ppmOI+B0offset);
                    negB0corr_ppm = polyval(npfit, ppmOI-B0offset);
                    
                    pimgout(i,j) = posB0corr_ppm;
                    nimgout(i,j) = negB0corr_ppm;
                    
                end
            end
        end
        
    case 2 % NW fast version of new, no more loops - fits entire version of matrix always
        pfit = NWpolyfitim(poly_deg,ppmlist,posimages);
        nfit = NWpolyfitim(poly_deg,ppmlist,negimages);
        
        pimgout = NWpolyvalim(pfit,ppmOI+B0map);
        nimgout = NWpolyvalim(nfit,ppmOI-B0map);
        
        % matches output to previous opts
        % not sure if better or worse though
        % comment out if desired
        pimgout = mask.*pimgout; 
        nimgout = mask.*nimgout;
        
end

end