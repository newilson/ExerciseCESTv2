function out = reprocessExerciseCESTwithROIs(in,filtstrength,reconType)
%
% out = reprocessExerciseCESTwithROIs(in)
%
% goes through ExerciseCEST processing using existing ROIs, allowing for
% non ROI-based changes

if nargin<3 || isempty(reconType)
    reconType = 'dicom'; % reconType can be 'dicom' or 'mat'
end

out = in; % initialize output to input

% UPDATE SETTINGS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    filtstrength = 'strong';
end
switch filtstrength
    case 'strong'
        out.dofilt = [128,96,32];
    case 'moderate'
        out.dofilt = [128,112,24];
    case 'weak'
        out.dofilt = [128,128,16];
    otherwise
        out.dofilt = false;
end

use_bounds = true;
lb1 = [0 20 4]; % increase, recovery, baseline
ub1 = [20 600 12];
lb2 = [0 20 4 4]; % increase, recovery, baseline, baseline
ub2 = [20 600 12 12];

if length(in.params.par0)==4
    lb1 = lb2; ub1 = ub2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Reference
for ii={'pre','post'}
    if strcmp(reconType,'mat')
        if strcmp(ii,'pre')
            pathname1.hdr = out.pre.refhdr;
            pathname1.refimage = out.pre.refimage;
            pathname1.dicomhdr = out.pre.refdicomhdr;
        else
            pathname1.hdr = out.post.refhdr;
            pathname1.refimage = out.post.refimage;
            pathname1.dicomhdr = out.post.refdicomhdr;
        end
    else % 'dicom'
        if strcmp(ii,'pre')
            pathname1 = fullfile(out.mainDir,out.subdirs.Refpre);
        else
            pathname1 = fullfile(out.mainDir,out.subdirs.Refpost);
        end
    end
    [refimage,hdr,dicomhdr,refimageb,~,biasscale] = readref(pathname1,false); % do not filter reference images
    if size(refimage,3)>1
        refimage = squeeze(refimage(:,:,1));
        refimageb = squeeze(refimageb(:,:,1));
    end
    if strcmp(ii,'pre')
        out.pre.refdicomhdr = dicomhdr;
        out.pre.refhdr = hdr;
        out.pre.refimage = refimage;
        out.pre.bias = biasscale;
        out.pre.refimageb = refimageb;
    else
        out.post.refdicomhdr = dicomhdr;
        out.post.refhdr = hdr;
        out.post.refimage = refimage;
        out.post.bias = biasscale;
        out.post.refimageb = refimageb;
    end
end

% Get pre v post registration transform
if out.par.isregister
    refpre = out.pre.refimageb.*out.pre.roimask;
    refpost = out.post.refimageb.*out.post.roimask;
    tform1 = imregcorr(refpre,refpost,'rigid');
    
    [optimizer,metric] = imregconfig('monomodal');
    tform = imregtform(refpre,refpost,'rigid',optimizer,metric,'InitialTransformation',tform1); % uses affine transform between pre/post exercise
    
    out.pre.tform = tform;
    
    % Apply transformation to ALL pre-exercise scans
    out.txtDisp.Value = 'Applying tranformations';
    out.pre.refimage = imwarp(out.pre.refimage,tform,'OutputView',imref2d(size(out.post.refimage)));
    out.pre.refimageb = imwarp(out.pre.refimageb,tform,'OutputView',imref2d(size(out.post.refimageb)));
    out.pre.roimask = out.post.roimask;
end

% B0
if out.par.iswassr
    for ii={'pre','post'}
        if strcmp(reconType,'mat')
            if strcmp(ii,'pre')
                pathname1.hdr = out.pre.WASSRhdr;
                images = [];
                for zz=1:size(out.pre.WASSRposimages,ndims(out.pre.WASSRposimages))
                    images = cat(3,images,out.pre.WASSRposimages(:,:,zz)); % 2d only right now
                    images = cat(3,images,out.pre.WASSRnegimages(:,:,zz));
                end
                pathname1.images = images;
                roimask = [];
                clear images
            else
                pathname1.hdr = out.post.WASSRhdr;
                images = [];
                for zz=1:size(out.post.WASSRposimages,ndims(out.post.WASSRposimages))
                    images = cat(3,images,out.post.WASSRposimages(:,:,zz)); % 2d only right now
                    images = cat(3,images,out.post.WASSRnegimages(:,:,zz));
                end
                pathname1.images = images;
                roimask = [];
                clear images
            end
        else % 'dicom'
            if strcmp(ii,'pre')
                pathname1 = fullfile(out.mainDir,out.subdirs.WASSRpre);
                roimask = out.pre.roimask;
            else
                pathname1 = fullfile(out.mainDir,out.subdirs.WASSRpost);
                roimask = out.post.roimask;
            end
        end
        
        if strcmp(pathname1,out.mainDir)
            B0map = zeros(size(roimask)); hdr = []; pars = []; posimages = []; negimages = [];
        else
            [B0map,hdr,pars,posimages,negimages,roimask] = readwassr(pathname1,1e6,roimask,out.dofilt); % arbitrary large number
        end
        
        if ~isequal(size(roimask),size(B0map))
            B0map = imresize(B0map,size(roimask));
            warndlg('B0map was resized')
        end
        
        if out.par.isregister
            B0map = imwarp(B0map,out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
        end
        
        roimask( abs(B0map) > out.params.b0tol ) = 0;
        roimask( B0map==0 ) = 0;
        
        if strcmp(ii,'pre')
            out.pre.B0map = B0map;
            out.pre.WASSRposimages = posimages;
            out.pre.WASSRnegimages = negimages;
            out.pre.WASSRhdr = hdr;
            if ~isempty(pars)
                out.pre.WASSRppmlist = pars.ppmlist;
                out.pre.nWASSRppm = length( pars.ppmlist);
                out.pre.WASSRb1 = pars.b1;
                out.pre.WASSRpw = pars.pw;
                out.pre.WASSRpw1 = pars.pw1;
                out.pre.WASSRdc = pars.dc;
            end
            out.pre.roimask = roimask;
        else
            out.post.B0map = B0map;
            out.post.WASSRposimages = posimages;
            out.post.WASSRnegimages = negimages;
            out.post.WASSRhdr = hdr;
            if ~isempty(pars)
                out.post.WASSRppmlist = pars.ppmlist;
                out.post.nWASSRppm = length( pars.ppmlist);
                out.post.WASSRb1 = pars.b1;
                out.post.WASSRpw = pars.pw;
                out.post.WASSRpw1 = pars.pw1;
                out.post.WASSRdc = pars.dc;
            end
            out.post.roimask = roimask;
        end
    end
else % TODO: add reconType 'mat' later
    for ii={'pre','post'}
        if strcmp(ii,'pre')
            pathname1 = fullfile(out.mainDir,out.subdirs.B0magpre);
            pathname2 = fullfile(out.mainDir,out.subdirs.B0phpre);
            roimask = out.pre.roimask;
        else
            pathname1 = fullfile(out.mainDir,out.subdirs.B0magpost);
            pathname2 = fullfile(out.mainDir,out.subdirs.B0phpost);
            roimask = out.post.roimask;
        end
        
        if strcmp(pathname2,out.mainDir)
            B0map = zeros(size(roimask)); hdr = []; ampl = []; phases = [];
        else
            [B0map, hdr, ampl, phases] = readB0gre(pathname1,pathname2,out.dofilt);
        end
        
        roimask( abs(B0map) > out.params.b0tol ) = 0;
        roimask( B0map==0 ) = 0;
        
        if ~isequal(size(roimask),size(B0map))
            B0map = imresize(B0map,size(roimask));
            warndlg('B0map was resized')
        end
        
        if out.par.isregister
            B0map = imwarp(B0map,out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
            for jj=1:size(ampl,3)
                ampl(:,:,jj) = imwarp(ampl(:,:,jj),out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
                phases(:,:,jj) = imwarp(phases(:,:,jj),out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
            end
        end
        
        if strcmp(ii,'pre')
            out.pre.B0GREA = ampl;
            out.pre.B0GREP = phases;
            if isempty(hdr)
                out.pre.B0GREhdr = [];
            else
                out.pre.B0GREhdr = hdr(1);
            end
            out.pre.B0map = B0map;
            out.pre.roimask = roimask;
        else
            out.post.B0GREA = ampl;
            out.post.B0GREP = phases;
            if isempty(hdr)
                out.post.B0GREhdr = [];
            else
                out.post.B0GREhdr = hdr(1);
            end
            out.post.B0map = B0map;
            out.post.roimask = roimask;
        end
    end
end

% B1
for ii={'pre','post'}
    if strcmp(reconType,'mat')
        if strcmp(ii,'pre')
            pathname1.hdr = out.pre.B1hdr;
            pathname1.images = out.pre.B1images;
            roimask = out.pre.roimask;
        else
            pathname1.hdr = out.post.B1hdr;
            pathname1.images = out.post.B1images;
            roimask = out.post.roimask;
        end
    else
        if strcmp(ii,'pre')
            pathname1 = fullfile(out.mainDir,out.subdirs.B1pre);
            roimask = out.pre.roimask;
        else
            pathname1 = fullfile(out.mainDir,out.subdirs.B1post);
            roimask = out.post.roimask;
        end
    end
    
    if strcmp(pathname1,out.mainDir)
        B1map = ones(size(roimask)); hdr = []; alpha = []; images = [];
    else
        [B1map,hdr,alpha,images,pathname1] = readB1(pathname1,out.dofilt);
    end
    
    if ~isequal(size(roimask),size(B1map))
        B1map = imresize(B1map,size(roimask));
        warndlg('B1map was resized')
    end
    
    roimask( B1map<1-out.params.b1tol ) = 0;
    roimask( B1map>1+out.params.b1tol ) = 0;
    
    if strcmp(ii,'pre')
        out.pre.B1map = B1map;
        out.pre.B1images = images;
        out.pre.alpha = alpha;
        out.pre.roimask = roimask;
        if isempty(hdr)
            out.pre.B1hdr = [];
        else
            out.pre.B1hdr = hdr(1);
        end
    else
        out.post.B1map = B1map;
        out.post.B1images = images;
        out.post.alpha = alpha;
        out.post.roimask = roimask;
        if isempty(hdr)
            out.post.B1hdr = [];
        else
            out.post.B1hdr = hdr(1);
        end
    end
end

% CEST images
for ii={'pre','post'}
    if strcmp(reconType,'mat')
        if strcmp(ii,'pre')
            pathname1.hdr = out.pre.CESThdr;
            images = [];
            for zz=1:size(out.pre.CESTposimages,ndims(out.pre.CESTposimages))
                images = cat(3,images,out.pre.CESTposimages(:,:,zz)); % 2d only right now
                images = cat(3,images,out.pre.CESTnegimages(:,:,zz));
            end
            pathname1.images = images;
            clear images
        else
            pathname1.hdr = out.post.CESThdr;
            images = [];
            for zz=1:size(out.post.CESTposimages,ndims(out.post.CESTposimages))
                images = cat(3,images,out.post.CESTposimages(:,:,zz)); % 2d only right now
                images = cat(3,images,out.post.CESTnegimages(:,:,zz));
            end
            pathname1.images = images;
            clear images
        end
    else % 'dicom'
        if strcmp(ii,'pre')
            if isfield(out.params,'multiscans') && out.params.multiscans
                pathname1 = out.subdirs.CESTprepaths;
            else
                pathname1 = fullfile(out.mainDir,out.subdirs.CESTpre);
            end
        else
            if isfield(out.params,'multiscans') && out.params.multiscans
                pathname1 = out.subdirs.CESTpostpaths;
            else
                pathname1 = fullfile(out.mainDir,out.subdirs.CESTpost);
            end
        end
    end
    
    if strcmp(pathname1,out.mainDir)
        posimages = []; negimages = []; hdr = []; pars = [];
    else
        [posimages, negimages, hdr, pars] = readCEST(pathname1,out.dofilt);
    end
    
    if out.par.isregister
        for jj=1:size(posimages,3)
            posimages(:,:,jj) = imwarp(posimages(:,:,jj),out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
            negimages(:,:,jj) = imwarp(negimages(:,:,jj),out.pre.tform,'OutputView',imref2d(size(out.post.refimage)));
        end
    end
    
    if strcmp(ii,'pre')
        out.pre.CESTpw = pars.CESTpw;
        out.pre.CESTb1 = pars.CESTb1;
        out.pre.CESTpw1 = pars.CESTpw1;
        out.pre.CESTdc = pars.CESTdc;
        out.pre.reqreps = pars.reqreps;
        out.pre.mocoreps = pars.mocoreps;
        out.pre.CESTppmlist = pars.CESTppmlist;
        out.pre.CESTposimages = posimages;
        out.pre.CESTnegimages = negimages;
        out.pre.CESThdr = hdr(1);
    else
        out.post.CESTpw = pars.CESTpw;
        out.post.CESTb1 = pars.CESTb1;
        out.post.CESTpw1 = pars.CESTpw1;
        out.post.CESTdc = pars.CESTdc;
        out.post.reqreps = pars.reqreps;
        out.post.mocoreps = pars.mocoreps;
        out.post.CESTppmlist = pars.CESTppmlist;
        out.post.CESTposimages = posimages;
        out.post.CESTnegimages = negimages;
        out.post.CESThdr = hdr(1);
    end
end

% Parameter Checking
if max(abs(out.pre.CESTppmlist-out.params.cestppm))>out.params.b0tol+1e-6 || max(abs(out.post.CESTppmlist-out.params.cestppm))>out.params.b0tol+1e-6
    warndlg('B0 tolerance is too large for the acquired offsets.')
end

% Corrected CEST
for ii={'pre','post'}
    if strcmp(ii,'pre')
        roimask = out.pre.roimask;
        posimages = out.pre.CESTposimages;
        negimages = out.pre.CESTnegimages;
        ppmlist = out.pre.CESTppmlist;
        reqreps = out.pre.reqreps;
        mocoreps = out.pre.mocoreps;
        B1map = out.pre.B1map;
        B0map = out.pre.B0map;
        CESTb1 = out.pre.CESTb1;
        refimage = out.pre.refimage;
    else
        roimask = out.post.roimask;
        posimages = out.post.CESTposimages;
        negimages = out.post.CESTnegimages;
        ppmlist = out.post.CESTppmlist;
        reqreps = out.post.reqreps;
        mocoreps = out.post.mocoreps;
        B1map = out.post.B1map;
        B0map = out.post.B0map;
        CESTb1 = out.post.CESTb1;
        refimage = out.post.refimage;
    end
    
    nCEST = length(ppmlist);
    rawCESTmap = zeros(size(posimages,1),size(posimages,2),reqreps*mocoreps);
    imageP = 0*rawCESTmap; imageN = 0*rawCESTmap;
    imageindex = find( abs(ppmlist - out.params.cestppm) == min(abs(ppmlist - out.params.cestppm)) );
    for cc=1:reqreps*mocoreps
        image1 = negimages(:,:,imageindex+(cc-1)*nCEST);
        image2 = posimages(:,:,imageindex+(cc-1)*nCEST);
        if out.par.norm==0 % Neg
            normppm = -out.params.cestppm;
            normimage = image1;
        elseif out.par.norm==2 % MT
            [fname,pname] = uigetfile({[out.mainDir filesep '*.dcm'];[out.mainDir filesep '*.ima']},...
                'Choose dicom for normalization');
            [hdr,normimage] = readdicomfiles2d(fullfile(pname,fname));
            normppm = inputdlg('Enter saturation offset in ppm ','MT offset');
        else % S0
            normppm = NaN;
            normimage = refimage;
        end
        image3 = 100*(image1-image2) ./ (normimage+eps) .* roimask;
        if out.fit.isfilt
            image3 = anisodiff(image3,20, 50, 0.03, 1);
        end
        rawCESTmap(:,:,cc) = image3;
        imageP(:,:,cc) = image2;
        imageN(:,:,cc) = image1;
    end
    if strcmp(ii,'pre')
        out.pre.rawCESTmap = rawCESTmap;
        out.pre.normimage = normimage;
        out.pre.imageP = imageP;
        out.pre.imageN = imageN;
    else
        out.post.rawCESTmap = rawCESTmap;
        out.post.normimage = normimage;
        out.post.imageP = imageP;
        out.post.imageN = imageN;
    end
        
    b0corrmaps = zeros(size(rawCESTmap));
    b0corrp = zeros(size(rawCESTmap));
    b0corrn = zeros(size(rawCESTmap));
    b0b1corrmaps = zeros(size(rawCESTmap));
    for cc = 1:reqreps*mocoreps
        pimg = posimages(:,:,1+(cc-1)*nCEST:cc*nCEST);
        nimg = negimages(:,:,1+(cc-1)*nCEST:cc*nCEST);
        if ( min(ppmlist) > 0.01 )
            B0mode = 'matlab';
            switch B0mode
                case 'mex'
                    [pout,nout] = mexB0correctedCESTimglfit(out.params.cestppm,ppmlist,pimg,nimg,B0map,double(roimask));
                case 'matlab'
                    poly_deg = 2;
                    [pout,nout] = B0correctedCEST2Dimg_ASNW(out.params.cestppm,ppmlist,pimg,nimg,B0map,roimask,poly_deg);
            end
        else
            nppm = length(ppmlist);
            [ppmsort,indexsort] = sort(ppmlist);
            pimgin = pimg*0;
            nimgin = nimg*0;
            for aa = 1:nppm
                pimgin(:,:,aa) = pimg(:,:,indexsort(aa));
                nimgin(:,:,aa) = nimg(:,:,indexsort(aa));
            end
            [pout,nout] = mexB0corr_interp_map_zspecinteg(out.params.cestppm,ppmsort,pimgin,nimgin,B0map,double(roimask));
        end
        if out.fit.isfilt
            B0corrp = anisodiff(pout,20,50,0.03,1);
            B0corrn = anisodiff(nout,20,50,0.03,1);
        else
            B0corrp = pout;
            B0corrn = nout;
        end
        if out.par.norm==0 % Neg
            normimage = B0corrn; % update to use the fitted negative offset
        end
        
        B0corrmap = 100 * (B0corrn - B0corrp)./(normimage+eps);
        B0corrmap(~isfinite(B0corrmap)) = 0;
        b0corrmaps(:,:,cc) = B0corrmap;
        b0corrp(:,:,cc) = B0corrp;
        b0corrn(:,:,cc) = B0corrn;
        
        
        B1corrtype = 0;
        switch B1corrtype
            case 0 % inverse
                B1corrmap = roimask .* (B0corrmap ./ (B1map+eps));
                B1corrmap(~isfinite(B1corrmap)) = 0;
            case 1 % calibration data
                if (handles.CF > 250.0)
                    B1calib = [17.54 -2.36]; % 7T
                else
                    B1calib = [ 8.29 -1.98]; % 3T
                end
                b1scale = CESTb1 * 0.612/42.5756;
                B1CR = B1map .* roimask * b1scale; % Abolute B1rms
                CESTmax = b1scale * (B1calib(2)*b1scale + B1calib(1));
                B1corrmap = B0corrmap + ( CESTmax - (B1calib(2)*(B1CR.^2) + B1calib(1)*B1CR) );
        end
        if out.fit.isfilt
            B1corrmap = anisodiff(B1corrmap,20,50,0.03,1);
        end
        b0b1corrmaps(:,:,cc) = B1corrmap;
    end
    if strcmp(ii,'pre')
        out.pre.B0corrCEST = b0corrmaps;
        out.pre.B0B1corrCEST = b0b1corrmaps;
        out.pre.B0corrp = b0corrp;
        out.pre.B0corrn = b0corrn;
    else
        out.post.B0corrCEST = b0corrmaps;
        out.post.B0B1corrCEST = b0b1corrmaps;
        out.post.B0corrp = b0corrp;
        out.post.B0corrn = b0corrn;
    end
    
    
end
out.pre.CESTmap = out.pre.B0B1corrCEST;
out.post.CESTmap = out.post.B0B1corrCEST;

% Fit recovery curves
refpre = out.pre.refimageb .* out.pre.roimask;
refpost = out.post.refimageb .* out.post.roimask;

if ~isfield(out.params,'myfun')
    myfun = [];
    par0 = [];
    parnames = [];
    skip_fit = true;
    out.fit.skip_fit = true;
else
    myfun = out.params.myfun;
    par0 = out.params.par0;
    skip_fit = false;
    out.fit.skip_fit = false;
    
    if isfield(out.params,'parnames')
        parnames = out.params.parnames;
    else
        parnames = [];
    end
    npars = length(par0);
    if isequal(npars,length(parnames))
        parnames{npars+1} = 'Voxels Post';
        parnames{npars+2} = 'Voxels Pre';
    end
    out.params.parnames = parnames;
end

cestpre = out.pre.CESTmap;
cestpost = out.post.CESTmap;
cest = cat(3,cestpre,cestpost);
npre = size(cestpre,3);
npost = size(cestpost,3);
nx = size(cestpre,1);
ny = size(cestpre,2);

roipre = out.pre.rois;
roipost = out.post.rois;
nrois = size(roipre,3);
timepoints = out.fit.timepoints;

nthr_abs = 0; nthr_per = 0;
roipostth_abs = []; roipostth_per = [];
if isfield(out.params,'roithresh')
    nthr_abs = length(out.params.roithresh);
    roipostth_abs = zeros([size(refpre), nrois, nthr_abs]);
    for kk=1:nthr_abs
        for ii = 1:nrois
            if ~isfield(out.params,'roithreshopt') || out.params.roithreshopt==0
                cestdiffpost = cest(:,:,npre+1)-cest(:,:,end); % post compared to recovery point-by-point
            else
                ind = find(roipost(:,:,ii));
                [ind1,ind2] = ind2sub([nx,ny],ind);
                cestdiffpost = cest(:,:,npre+1)-mean(col(cest(ind1,ind2,1:npre))); % post compared to mean pre ROI
            end
            roipostth_abs(:,:,ii,kk) = (cestdiffpost>out.params.roithresh(kk)) .* logical(roipost(:,:,ii));
        end
    end
end
if isfield(out.params,'roithreshpercent')
    nthr_per = length(out.params.roithreshpercent);
    roipostth_per = zeros([size(refpre), nrois, nthr_per]);
    for kk=1:nthr_per
        for ii=1:nrois
            if ~isfield(out.params,'roithreshpercentopt') || out.params.roithreshpercentopt==0
                cestdiffpost = 100*(cest(:,:,npre+1)-cest(:,:,end))./cest(:,:,end);
            else
                ind = find(roipost(:,:,ii));
                [ind1,ind2] = ind2sub([nx,ny],ind);
                cestdiffpost = 100*(cest(:,:,npre+1)-mean(col(cest(ind1,ind2,1:npre))))/mean(col(cest(ind1,ind2,1:npre)));
            end
            roipostth_per(:,:,ii,kk) = (cestdiffpost>out.params.roithreshpercent(kk)) .* logical(roipost(:,:,ii));
        end
    end
end
nthr = max(nthr_abs+nthr_per,1);
out.fit.nthr = nthr;
out.fit.nthr_abs = nthr_abs;
out.fit.nthr_per = nthr_per;
roipostth = cat(4,roipostth_abs,roipostth_per);

weighting = out.params.weights;
nweights = length(weighting);
means = ones(nrois,npre+npost,nthr); stds = ones(nrois,npre+npost,nthr);
a = zeros(length(par0),nrois,nweights,nthr); r = zeros(npre+npost,nrois,nweights,nthr); J = zeros(npre+npost,length(par0),nrois,nweights,nthr);
cov = zeros(length(par0),length(par0),nrois,nweights,nthr); mse = zeros(1,nrois,nweights,nthr);
fitdata = zeros(length(timepoints),nrois,nweights,nthr);
ci95 = zeros(length(par0),2,nrois,nweights,nthr); ci90 = ci95;
res = 0*r; 
resnorm = zeros(nrois,nweights,nthr);
out.fit.nfit = 1;
if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
    out.fit.nfit = 2;
    par0_2 = out.params.par0_2;
    myfun2 = out.params.myfun_2;
    parnames2 = out.params.parnames_2;
    a2 = zeros(length(par0_2),nrois,nweights,nthr); r2 = zeros(npre+npost,nrois,nweights,nthr); J2 = zeros(npre+npost,length(par0_2),nrois,nweights,nthr);
    cov2 = zeros(length(par0_2),length(par0_2),nrois,nweights,nthr); mse2 = zeros(1,nrois,nweights,nthr);
    fitdata2 = zeros(length(timepoints),nrois,nweights,nthr);
    ci95_2 = zeros(length(par0_2),2,nrois,nweights,nthr); ci90_2 = ci95_2;
    res2 = 0*r2;
    resnorm2 = zeros(nrois,nweights,nthr);
end
npts_interp = 100;
firstpost = zeros(nrois,nthr); Tt50sp = zeros(nrois,nthr); Tt50sp2 = zeros(nrois,nthr); baseline = zeros(nrois,1);
means_sp = zeros(nrois,npts_interp,nthr); means_sp2 = zeros(nrois,npts_interp,nthr);
warning('off','stats:nlinfit:IllConditionedJacobian');
for kk=1:nthr % number of thresholds
    for ii=1:nrois % number of rois
        for jj=1:npre+npost % number of timepoints
            cestpoint = cest(:,:,jj);
            if jj<npre+1
                cestroi = cestpoint(roipre(:,:,ii)>0);
            else
                cestroi = cestpoint(roipost(:,:,ii)>0);
                if isfield(out.params,'roithresh') || isfield(out.params,'roithreshpercent')
                    cestroith = cestpoint(roipostth(:,:,ii,kk)>0);
                else
                    cestroith = [];
                end
            end
            if jj>=npre+1 && ~isempty(cestroith)
                means(ii,jj,kk) = mean(cestroith);
                stds(ii,jj,kk) = std(cestroith);
            else % use whole roi for pre exercise
                means(ii,jj,kk) = mean(cestroi);
                stds(ii,jj,kk) = std(cestroi);
            end
        end
        if ~skip_fit
            for zz=1:nweights
                switch weighting{zz}%button
                    case '1'
                        weights = ones(size(means));
                    case '1/N'
                        weights = 1./stds;
                    case 'S'
                        weights = means;
                    case 'S/N'
                        weights = means./stds;
                    otherwise
                        weights = ones(size(means));
                end
                weights = abs(weights);
                for ff=1:2
                    try
                        if ff==1
                            if ~use_bounds
                                [a(:,ii,zz,kk),r(:,ii,zz,kk),J(:,:,ii,zz,kk),cov(:,:,ii,zz,kk),mse(:,ii,zz,kk)] = nlinfit(timepoints,means(ii,:,kk),myfun,par0,'weights',weights(ii,:,kk));
                                ci95(:,:,ii,zz,kk) = nlparci(a(:,ii,zz,kk),r(:,ii,zz,kk),'covar',cov(:,:,ii,zz,kk),'alpha',0.05);
                                ci90(:,:,ii,zz,kk) = nlparci(a(:,ii,zz,kk),r(:,ii,zz,kk),'covar',cov(:,:,ii,zz,kk),'alpha',0.1);
                            else
                                options = optimoptions(@lsqnonlin);
                                options.Display = 'off';
                                mycost = @(p) sqrt(weights(ii,:,kk)) .* (myfun(p,timepoints) - means(ii,:,kk));
                                [a(:,ii,zz,kk),resnorm(ii,zz,kk),res(:,ii,zz,kk)] = lsqnonlin(mycost,par0,lb1,ub1,options);
                            end
                            fitdata(:,ii,zz,kk) = myfun(a(:,ii,zz,kk),timepoints);
                        elseif ff==2
                            % NW temp
                            if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
                                if ~use_bounds
                                    [a2(:,ii,zz,kk),r2(:,ii,zz,kk),J2(:,:,ii,zz,kk),cov2(:,:,ii,zz,kk),mse2(:,ii,zz,kk)] = nlinfit(timepoints,means(ii,:,kk),myfun2,par0_2,'weights',weights(ii,:,kk));
                                    ci95_2(:,:,ii,zz,kk) = nlparci(a2(:,ii,zz,kk),r2(:,ii,zz,kk),'covar',cov2(:,:,ii,zz,kk),'alpha',0.05);
                                    ci90_2(:,:,ii,zz,kk) = nlparci(a2(:,ii,zz,kk),r2(:,ii,zz,kk),'covar',cov2(:,:,ii,zz,kk),'alpha',0.1);
                                else
                                    options = optimoptions(@lsqnonlin);
                                    options.Display = 'off';
                                    mycost2 = @(p) sqrt(weights(ii,:,kk)) .* (myfun2(p,timepoints) - means(ii,:,kk));
                                    [a2(:,ii,zz,kk),resnorm2(ii,zz,kk),res2(:,ii,zz,kk)] = lsqnonlin(mycost2,par0_2,lb2,ub2,options);
                                end
                                fitdata2(:,ii,zz,kk) = myfun2(a2(:,ii,zz,kk),timepoints);
                            end
                        end
                    catch errorObj
                        if ff==1
                            a(:,ii,zz,kk) = NaN;
                        elseif ff==2
                            a2(:,ii,zz,kk) = NaN;
                        end
                        disp(['fit error with ROI' num2str(ii) ', threshold ' num2str(kk) ', weight ' num2str(zz) ', function ' num2str(ff)]);
                    end
                end
            end
        end
        
        % Simple fit proxy variables
        means_vec = squeeze(means(ii,:,kk));
        timepoints_interp = linspace(timepoints(npre+1),timepoints(end),npts_interp);
        ft = fittype('smoothingspline');
        opts = fitoptions('Method','SmoothingSpline');
        opts.SmoothingParam = 1e-6;
        [xData,yData] = prepareCurveData(col(timepoints(npre+1:end)),col(means_vec(npre+1:end)));
        [fitresult,gof] = fit(xData,yData,ft,opts);
        means_sp(ii,:,kk) = feval(fitresult,timepoints_interp);
        means_sp2(ii,:,kk) = interp1(timepoints(npre+1:end),means_vec(npre+1:end),timepoints_interp,'spline');
        baseline(ii) = mean(means_vec(1:npre));
        firstpost(ii,kk) = means_vec(npre+1);
        halfrecov = (firstpost(ii,kk)-baseline(ii))/2 + baseline(ii);
        [~,indsp] = min(abs(squeeze(means_sp(ii,:,kk))-halfrecov));
        Tt50sp(ii,kk) = timepoints_interp(indsp);
        [~,indsp2] = min(abs(squeeze(means_sp2(ii,:,kk))-halfrecov));
        Tt50sp2(ii,kk) = timepoints_interp(indsp2);
    end
end
warning('on','stats:nlinfit:IllConditionedJacobian');

% Data-based slopes
slope2 = ( means(:,npre+2,:)-means(:,npre+1,:) ) ./ ( timepoints(npre+2)-timepoints(npre+1) ); % first 2 points only, no regression
slope_sp = ( means_sp(:,2,:) - means_sp(:,1,:) ) ./ ( timepoints_interp(2)-timepoints_interp(1) );

p3 = NWpolyfitim(1,timepoints(npre+1:npre+3),permute(means(:,npre+1:npre+3,:),[1 3 2]));
p4 = NWpolyfitim(1,timepoints(npre+1:npre+4),permute(means(:,npre+1:npre+4,:),[1 3 2]));
p5 = NWpolyfitim(1,timepoints(npre+1:npre+5),permute(means(:,npre+1:npre+5,:),[1 3 2]));
p10 = NWpolyfitim(1,timepoints_interp(1:10),permute(means_sp(:,1:10,:),[1 3 2]));
slope3 = squeeze(p3(:,:,1));
slope4 = squeeze(p4(:,:,1));
slope5 = squeeze(p5(:,:,1));
slope10 = squeeze(p10(:,:,1));

maxsig = means+stds;
maxsig = max(maxsig(:));
out.fit.nrois = nrois;
if ~skip_fit
    out.fit.pars = a;
    if ~use_bounds
        out.fit.res = r;
        out.fit.J = J;
        out.fit.cov = cov;
        out.fit.mse = mse;
        out.fit.confidence_intervals95 = ci95;
        out.fit.confidence_intervals90 = ci90;
    else
        out.fit.resnorm = resnorm;
        out.fit.res = res;
        out.fit.lb1 = lb1;
        out.fit.ub1 = ub1;
        out.fit.lb2 = lb2;
        out.fit.ub2 = ub2;
    end
    out.fit.fitdata = fitdata;
    out.fit.nweights = nweights;
    % NW temp
    if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
        out.params.parnames_2{length(out.params.par0_2)+1} = 'Voxels Post';
        out.params.parnames_2{length(out.params.par0_2)+2} = 'Voxels Pre';
        out.fit.pars_2 = a2;
        if ~use_bounds
            out.fit.res_2 = r2;
            out.fit.J_2 = J2;
            out.fit.cov_2 = cov2;
            out.fit.mse_2 = mse2;
            out.fit.confidence_intervals95_2 = ci95_2;
            out.fit.confidence_intervals90_2 = ci90_2;
        else
            out.fit.resnorm_2 = resnorm2;
            out.fit.res_2 = res2;
        end
        out.fit.fitdata_2 = fitdata2;
    end
end
out.fit.means = means;
out.fit.stds = stds;
out.fit.timepoints = timepoints;
out.fit.plotmax = maxsig;
out.fit.roipix.tot = sum(out.fit.roipix.pre);
if isfield(out.params,'roithresh')
    out.post.roisth = roipostth;
    for kk=1:nthr
        out.fit.roipix.totpost(kk) = sum(col(roipostth(:,:,:,kk)));
        for ii=1:nrois
            temp = roipostth(:,:,ii,kk);
            out.fit.roipix.post(kk,ii) = sum(temp(:));
        end
    end
end
if ~isfield(out.fit.roipix,'totpost'), out.fit.roipix.totpost = sum(out.fit.roipix.post); end
out.fit.roipix.totfrac = out.fit.roipix.totpost/out.fit.roipix.tot;
out.fit.nomodelpars.baseline = baseline;
out.fit.nomodelpars.firstpost = firstpost;
out.fit.nomodelpars.increase = firstpost-baseline;
out.fit.nomodelpars.Tt50sp = Tt50sp2;
out.fit.nomodelpars.slope2 = squeeze(slope2);
out.fit.nomodelpars.slope3 = slope3;
out.fit.nomodelpars.slope4 = slope4;
out.fit.nomodelpars.slope5 = slope5;
out.fit.nomdoelpars.means = means_sp2;
out.fit.sp.means = means_sp;
out.fit.sp.slope10 = slope10;
out.fit.sp.firstpoint = squeeze(means_sp(:,1,:));
out.fit.sp.slope2 = squeeze(slope_sp);
out.fit.sp.baseline = baseline;
out.fit.sp.increase = out.fit.sp.firstpoint-out.fit.sp.baseline;
out.fit.sp.Tt50 = Tt50sp;



% Global ROI if applicable
if out.params.globalroi
    gloroipre = out.pre.gloroi;
    gloroipost = out.post.gloroi;
    if isfield(out.params,'roithresh')
        roipostth_abs = zeros([size(out.post.roimask),nthr_abs]);
        for kk=1:nthr_abs
            if ~isfield(out.params,'roithreshopt') || out.params.roithreshopt==0
                cestdiffpost = cest(:,:,npre+1)-cest(:,:,end); % post compared to recovery
            else
                ind = find(gloroipost);
                [ind1,ind2] = ind2sub([nx,ny],ind);
                cestdiffpost = cest(:,:,npre+1)-mean(col(cest(ind1,ind2,1:npre))); % post compared to mean pre
            end
            roipostth_abs(:,:,kk) = logical(cestdiffpost>out.params.roithresh(kk)) .* logical(gloroipost);
        end
    end
    if isfield(out.params,'roithreshpercent')
        roipostth_per = zeros([size(out.post.roimask),nthr_per]);
        for kk=1:nthr_per
            if ~isfield(out.params,'roithreshpercentopt') || out.params.roithreshpercentopt==0
                cestdiffpost = 100*(cest(:,:,npre+1)-cest(:,:,end))./cest(:,:,end); % post compared to recovery
            else
                ind = find(gloroipost);
                [ind1,ind2] = ind2sub([nx,ny],ind);
                cestdiffpost = 100*(cest(:,:,npre+1)-mean(col(cest(ind1,ind2,1:npre))))/mean(col(cest(ind1,ind2,1:npre))); % post compared to mean pre
            end
            roipostth_per(:,:,kk) = logical(cestdiffpost>out.params.roithreshpercent(kk)) .* logical(gloroipost);
        end
    end
    roipostth = cat(3,roipostth_abs,roipostth_per);
    
    means = ones(npre+npost,nthr); stds = ones(npre+npost,nthr);
    a = zeros(length(par0),nweights,nthr); r = zeros(npre+npost,nweights,nthr); J = zeros(npre+npost,length(par0),nweights,nthr);
    cov = zeros(length(par0),length(par0),nweights,nthr); mse = zeros(1,nweights,nthr);
    fitdata = zeros(length(timepoints),nweights,nthr);
    ci95 = zeros(length(par0),2,nweights,nthr); ci90 = ci95;
    res = 0*r; 
    resnorm = zeros(nweights,nthr);
    if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
        a2 = zeros(length(par0_2),nweights,nthr); r2 = zeros(npre+npost,nweights,nthr); J2 = zeros(npre+npost,length(par0_2),nweights,nthr);
        cov2 = zeros(length(par0_2),length(par0_2),nweights,nthr); mse2 = zeros(1,nweights,nthr);
        fitdata2 = zeros(length(timepoints),nweights,nthr);
        ci95_2 = zeros(length(par0_2),2,nweights,nthr); ci90_2 = ci95_2;
        res2 = 0*r2; 
        resnorm2 = zeros(nweights,nthr);
    end
    firstpost = zeros(1,nthr); Tt50sp = zeros(1,nthr); Tt50sp2 = zeros(1,nthr);
    means_sp = zeros(npts_interp,nthr); means_sp2 = zeros(npts_interp,nthr);
    warning('off','stats:nlinfit:IllConditionedJacobian');
    for kk=1:nthr
        for jj=1:npre+npost
            cestpoint = cest(:,:,jj);
            if jj<npre+1
                cestroi = cestpoint(gloroipre>0);
            else
                cestroi = cestpoint(gloroipost>0);
                if isfield(out.params,'roithresh') || isfield(out.params,'roithreshpercent')
                    cestroith = cestpoint(roipostth(:,:,kk)>0);
                else
                    cestroith = [];
                end
            end
            if jj>=npre+1 && ~isempty(cestroith)
                means(jj,kk) = mean(cestroith);
                stds(jj,kk) = std(cestroith);
            else % use whole roi for pre exercise
                means(jj,kk) = mean(cestroi);
                stds(jj,kk) = std(cestroi);
            end
        end
        if ~skip_fit
            for zz=1:nweights
                switch zz%button
                    case 1%'1'
                        weights = ones(size(means));
                    case 2%'1/N'
                        weights = 1./stds;
                    case 3%'S'
                        weights = means;
                    case 4%'S/N'
                        weights = means./stds;
                end
                weights = abs(weights);
                for ff=1:2
                    try
                        if ff==1
                            if ~use_bounds
                                [a(:,zz,kk),r(:,zz,kk),J(:,:,zz,kk),cov(:,:,zz,kk),mse(:,zz,kk)] = nlinfit(timepoints,means(:,kk).',myfun,par0,'weights',weights(:,kk).');
                                ci95(:,:,zz,kk) = nlparci(a(:,zz,kk),r(:,zz,kk),'covar',cov(:,:,zz,kk),'alpha',0.05);
                                ci90(:,:,zz,kk) = nlparci(a(:,zz,kk),r(:,zz,kk),'covar',cov(:,:,zz,kk),'alpha',0.1);
                            else
                                options = optimoptions(@lsqnonlin);
                                options.Display = 'off';
                                mycost = @(p) sqrt(weights(:,kk))' .* (myfun(p,timepoints) - means(:,kk)');
                                [a(:,zz,kk),resnorm(zz,kk),res(:,zz,kk)] = lsqnonlin(mycost,par0,lb1,ub1,options);
                            end
                            fitdata(:,zz,kk) = myfun(a(:,zz,kk),timepoints);
                        elseif ff==2% NW temp
                            if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
                                if ~use_bounds
                                    [a2(:,zz,kk),r2(:,zz,kk),J2(:,:,zz,kk),cov2(:,:,zz,kk),mse2(:,zz,kk)] = nlinfit(timepoints,means(:,kk).',myfun2,par0_2,'weights',weights(:,kk).');
                                    ci95_2(:,:,zz,kk) = nlparci(a2(:,zz,kk),r2(:,zz,kk),'covar',cov2(:,:,zz,kk),'alpha',0.05);
                                    ci90_2(:,:,zz,kk) = nlparci(a2(:,zz,kk),r2(:,zz,kk),'covar',cov2(:,:,zz,kk),'alpha',0.1);
                                else
                                    options = optimoptions(@lsqnonlin);
                                    options.Display = 'off';
                                    mycost2 = @(p) sqrt(weights(:,kk))' .* (myfun2(p,timepoints) - means(:,kk)');
                                    [a2(:,zz,kk),resnorm2(zz,kk),res2(:,zz,kk)] = lsqnonlin(mycost2,par0_2,lb2,ub2,options);
                                end
                                fitdata2(:,zz,kk) = myfun2(a2(:,zz,kk),timepoints);
                            end
                        end
                    catch errorObj
                        if ff==1
                            a(:,zz,kk) = NaN;
                        elseif ff==2
                            a2(:,zz,kk) = NaN;
                        end
                        if length(findall(0,'Type','figure','Name','Error'))==0
                            errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
                        end
                        
                        disp(['fit error with Global, threshold ' num2str(kk) ', weight ' num2str(zz) ', function ' num2str(ff)]);
                        
                    end
                end
            end
        end
        
        % Simple fit proxy variables
        means_vec = squeeze(means(:,kk));
        timepoints_interp = linspace(timepoints(npre+1),timepoints(end),100);
        ft = fittype('smoothingspline');
        opts = fitoptions('Method','SmoothingSpline');
        opts.SmoothingParam = 1e-6;
        [xData,yData] = prepareCurveData(col(timepoints(npre+1:end)),col(means_vec(npre+1:end)));
        [fitresult,gof] = fit(xData,yData,ft,opts);
        means_sp(:,kk) = feval(fitresult,timepoints_interp);
        means_sp2(:,kk) = interp1(timepoints(npre+1:end),means_vec(npre+1:end),timepoints_interp,'spline');
        baseline = mean(means_vec(1:npre));
        firstpost(kk) = means_vec(npre+1);
        halfrecov = (firstpost(kk)-baseline)/2 + baseline;
        [~,indsp] = min(abs(means_sp(:,kk)-halfrecov));
        Tt50sp(kk) = timepoints_interp(indsp);
        [~,indsp2] = min(abs(means_sp2(:,kk)-halfrecov));
        Tt50sp2(kk) = timepoints_interp(indsp2);
   
    end
    warning('on','stats:nlinfit:IllConditionedJacobian');
    
    % Data-based slopes
    slope2 = ( means(npre+2,:)-means(npre+1,:) ) ./ ( timepoints(npre+2)-timepoints(npre+1) ); % first 2 points only, no regression
    slope_sp = ( means_sp(2,:) - means_sp(1,:) ) ./ ( timepoints_interp(2)-timepoints_interp(1) );

    p3 = NWpolyfitim(1,timepoints(npre+1:npre+3),permute(means(npre+1:npre+3,:),[2 3 1]));
    p4 = NWpolyfitim(1,timepoints(npre+1:npre+4),permute(means(npre+1:npre+4,:),[2 3 1]));
    p5 = NWpolyfitim(1,timepoints(npre+1:npre+5),permute(means(npre+1:npre+5,:),[2 3 1]));
    p10 = NWpolyfitim(1,timepoints_interp(1:10),permute(means_sp(1:10,:),[2 3 1]));
    slope3 = squeeze(p3(:,1));
    slope4 = squeeze(p4(:,1));
    slope5 = squeeze(p5(:,1));
    slope10 = squeeze(p10(:,1));
    
    if ~skip_fit
        out.fit.glo.pars = a;
        if ~use_bounds
            out.fit.glo.res = r;
            out.fit.glo.J = J;
            out.fit.glo.cov = cov;
            out.fit.glo.mse = mse;
            out.fit.glo.confidence_intervals95 = ci95;
            out.fit.glo.confidence_intervals90 = ci90;
        else
            out.fit.glo.resnorm = resnorm;
            out.fit.glo.res = res;
        end
        out.fit.glo.fitdata = fitdata;
        % NW temp
        if isfield(out.params,'myfun_2') && isfield(out.params,'par0_2') && isfield(out.params,'parnames_2')
            out.fit.glo.pars_2 = a2;
            if ~use_bounds
                out.fit.glo.res_2 = r2;
                out.fit.glo.J_2 = J2;
                out.fit.glo.cov_2 = cov2;
                out.fit.glo.mse_2 = mse2;
                out.fit.glo.confidence_intervals95_2 = ci95_2;
                out.fit.glo.confidence_intervals90_2 = ci90_2;
            else
                out.fit.glo.resnorm_2 = resnorm2;
                out.fit.glo.res_2 = res2;
            end
            out.fit.glo.fitdata_2 = fitdata2;
        end
    end
    out.fit.glo.means = means;
    out.fit.glo.stds = stds;
    if isfield(out.params,'roithresh')
        out.post.gloroith = roipostth;
        for kk=1:nthr
            temp = roipostth(:,:,kk);
            out.fit.roipix.postglo(kk) = sum(temp(:));
        end
    end
    out.fit.glo.nomodelpars.baseline = baseline;
    out.fit.glo.nomodelpars.firstpost = firstpost;
    out.fit.glo.nomodelpars.increase = firstpost-baseline;
    out.fit.glo.nomodelpars.Tt50.spline = Tt50sp2;
    out.fit.glo.nomodelpars.slope2 = squeeze(slope2);
    out.fit.glo.nomodelpars.slope3 = slope3;
    out.fit.glo.nomodelpars.slope4 = slope4;
    out.fit.glo.nomodelpars.slope5 = slope5;    
    out.fit.glo.nomodelpars.means = means_sp2;
    out.fit.glo.sp.means = means_sp;
    out.fit.glo.sp.slope10 = slope10;
    out.fit.glo.sp.firstpoint = squeeze(means_sp(:,1,:));
    out.fit.glo.sp.slope2 = squeeze(slope_sp);
    out.fit.glo.sp.baseline = baseline;
    out.fit.glo.sp.increase = out.fit.glo.sp.firstpoint-baseline;
    out.fit.glo.sp.Tt50 = Tt50sp;


end