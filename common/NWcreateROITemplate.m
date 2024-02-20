function [pos,mask] = NWcreateROITemplate(im,nrois,savename,roinames)
%
% im is the reference image - must be the same size as the acquired image
% nrois is the number of rois or a cell array of initial positions
% savename is the name of the saved .mat file

if nargin<2, nrois = 1; end
if nargin<3, 
    savename = [];
%     savename = 'ROITemplate'; 
end
if isnumeric(nrois)
    opt = 0;
elseif iscell(nrois)
    opt = 1;
    init_pos = nrois;
    nrois = length(init_pos);
else
    error('incorrect number of rois')
end
if nargin<4 || ~isequal(length(roinames),nrois)
    roinames = [];
end

hfig = figure;
ax = axes(hfig);
imagesc(ax,im);
colormap bone

mask = zeros([size(im) nrois]);
for ii=1:nrois
    if isempty(roinames)
        hbox = msgbox(['Choose ROI ' num2str(ii)]);
    else
        hbox = msgbox(['Choose ' roinames{ii}]);
    end
    uiwait(hbox)
    switch opt
        case 0
            h{ii} = impoly(ax);
        case 1
            h{ii} = impoly(ax,init_pos{ii});
    end
    setColor(h{ii},'r')
    pos{ii} = wait(h{ii});
    mask(:,:,ii) = createMask(h{ii});
end
close(hfig)

if strfind(savename,'mat')
    savename = savename(1:end-4);
end

if ~isempty(savename)
    save(savename,'pos')
end