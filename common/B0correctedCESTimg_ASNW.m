%% B0correctedCEST2Dimg_AS.m functing created by @Anup, July 2017 
% New function for correcting CEST weighted images for B0 inhomogenity
% 'ppm' at which outpot image is required
%  ppmlist-- list of ppm values at which input posImgs  were acquired, 
%  posImg, images at ppmlist
%  NegImg, images at -ppmlist
%  B0 map,  mask
%  poly_deg  --defalt = 2, (it can be increased to 3 or 4 depending upon
%  nature of posimg intensit vs ppmlist curve shape.
% NW - modified for 3d 11/18/2024

function [pimgout,nimgout,pfit,nfit] = B0correctedCESTimg_ASNW(ppmOI,ppmlist,posimages,negimages,B0map,mask,poly_deg,opt)

if nargin<8 || isempty(opt), opt = 2; end

if opt~=2
    warning('only OPTION 2 supported')
    opt = 2;
end

ppmlist = squeeze(ppmlist)';

% NW - removed deprecated options (see B0correctedCEST2dimg_ASNW.m for
% original
switch 2        
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