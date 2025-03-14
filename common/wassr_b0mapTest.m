function [B0map, Zspec_int, ppm_int] = wassr_b0map(ppmlist, posimage, negimage, mask, b0ppmstep)
    % Inputs:
    % ppmlist    - List of ppm values (1D array)
    % posimage   - Positive image data (3D/4D array: H x W x [Z] x nppm)
    % negimage   - Negative image data (3D/4D array: H x W x [Z] x nppm)
    % mask       - Mask for the image (2D/3D array: H x W x [Z])
    % b0ppmstep  - Step size for B0 map calculation (default: 0.005)

    % Default step size
    if nargin < 5
        b0ppmstep = 0.005;
    end

    if ~isequal(size(posimage), size(negimage))
        error('Unequal image sizes. Cannot get B0 map.');
    end

    % Get dimensions
    if ndims(posimage) == 3
        [H, W, nppm] = size(posimage);
        Z = 1;
    elseif ndims(posimage) == 4
        [H, W, Z, nppm] = size(posimage);
    else
        error('Unexpected size. Cannot get B0 map.');
    end
    nelem = H * W * Z;

    if nargin < 4 || isempty(mask)
        mask = squeeze(true(H, W, Z));
    else
        mask = logical(mask);
    end
    si = size(mask);
    if Z == 1, si(3) = 1; end
    if ~isequal(si,[H,W,Z])
        error('wrong mask');
    end

    % Sort ppm list in ascending order
    [ppmlist_asc, sort_idx] = sort(ppmlist);
    minppm = ppmlist_asc(1);
    maxppm = ppmlist_asc(end);
    stepppm = (maxppm - minppm) / (nppm - 1);

    % Interpolation factor
    intfactor = floor(stepppm / b0ppmstep);
    nzsp1 = 2 * nppm - 1;
    nzsp2 = intfactor * (nzsp1 - 1) + 1;
    stepppm2 = stepppm / intfactor;
    rangezsp2 = floor(0.2 / stepppm2);

    % Generate interpolated ppm coordinates
    ppm_int = linspace(-maxppm, maxppm, nzsp2);
    ppm_int = ppm_int(:);

    % Vectorize
    posimage = reshape(posimage, [], nppm);
    negimage = reshape(negimage, [], nppm);
    mask = mask(:);

    % Initialize outputs
    B0map = zeros(nelem, 1); % column vector first
    Zspec_int = zeros(nelem, nzsp2);

    % Precompute spline coefficients
    xzsp = [-flip(ppmlist_asc); ppmlist_asc(2:end)];

    for ielem = 1:nelem
        if mask(ielem)
            % Extract z-spectrum for the current pixel
            zspp1 = posimage(ielem, :); zspp1 = zspp1(:);
            zspn1 = negimage(ielem, :); zspn1 = zspn1(:);

            % Sort z-spectrum
            zsp1 = [flip(zspn1(sort_idx)); zspp1(sort_idx(2:end))];

            % Find minimum of z-spectrum
            [~, min_idx] = min(zsp1);
            zspxmin = xzsp(min_idx);

            % Spline interpolation
            y2 = spline(xzsp, zsp1);
        
            % Evaluate spline
            zsp2 = ppval(y2, ppm_int);

            % Reflected z-spectrum
            zsp2r = flip(zsp2);

            % Find B0 shift using maximum symmetry algorithm
            temp1 = -zspxmin;
            ftemp2 = inf;
            minSSEindex = -rangezsp2;
            for i1 = -rangezsp2:rangezsp2
                temp2 = temp1 + (i1 * stepppm2);
                ftemp1 = maxSymWASSR(temp2, xzsp, zsp1, ppm_int, zsp2r, stepppm2);
                if ftemp1 < ftemp2
                    minSSEindex = i1;
                    ftemp2 = ftemp1;
                end
            end

            % Store B0 map value
            B0map(ielem) = -(temp1 + (minSSEindex * stepppm2));

            % Store interpolated z-spectrum
            Zspec_int(ielem, :) = zsp2;
        end
    end

    % Reshape to matrix
    B0map = reshape(B0map, H, W, Z);
    if Z > 1
        Zspec_int = reshape(Zspec_int, H, W, Z, nzsp2);
    else
        Zspec_int = reshape(Zspec_int, H, W, nzsp2);
    end
    ppm_int = ppm_int(:);
end

function SSE = maxSymWASSR(xval, xsp1, zsp1, xsp2, zsp2r, xstep2)
    % Maximum symmetry WASSR algorithm
    cf = 2 * xval; % Compute the correction factor
    A1 = xsp1 + cf; % Shift xsp1 by the correction factor

    % Define the valid range for A1
    if cf<0
        A2 = xsp1(1); % Lower bound
        A3 = cf + xsp1(end); % Upper bound
    else
        A2 = cf + xsp1(1); % Lower bound
        A3 = xsp1(end); % Upper bound
    end

    % Find valid indices where A1 is within the range [A2, A3]
    valid1 = (A1 > A2) & (A1 < A3);

    % Compute indices for zsp2r
    indices = round((A1 - xsp2(1)) / xstep2) + 1;

    % Ensure indices are within bounds of zsp2r
    valid2 = indices >= 1 & indices <= length(zsp2r);

    valid3 = valid1 & valid2;
    % Compute SSE only for valid indices
    if any(valid3)
        % Extract valid zsp1 and zsp2r values
        zsp1_valid = zsp1(valid3);
        zsp2r_valid = zsp2r(indices(valid3));

        % Compute the sum of squared differences (SSE)
        SSE = sum((zsp1_valid(:) - zsp2r_valid(:)).^2);

%         SSE = sum((zsp1(valid3) - zsp2r(indices(valid3))).^2);
    else
        % If no valid indices are found, return a large value (inf)
        SSE = inf;
    end
end
% 
% function SSE = maxSymWASSR(xval, xsp1, zsp1, xsp2, zsp2r, xstep2)
%     % Maximum symmetry WASSR algorithm
%     cf = 2 * xval;
%     A1 = xsp1 + cf;
%     A2 = cf + xsp1(1);
%     A3 = cf + xsp1(end);
% 
%     valid = (A1 > A2) & (A1 < A3);
%     indices = round((A1(valid) - xsp2(1)) / xstep2) + 1;
%     valid_indices = indices >= 1 & indices <= length(zsp2r);
% 
%     if any(valid_indices)
%         SSE = sum((zsp1(valid(valid_indices)) - zsp2r(indices(valid_indices))).^2);
%     else
%         SSE = inf;
%     end
% end