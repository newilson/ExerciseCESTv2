function NWimoverlay(Im,alpha,ax,colors)
%
% Im is 3d matrix of base image and overlay(s)
% colors is a [nx3] matrix of rgb
%
% Example
% im1 = phantom(128);
% im2 = 0*im1;
% im2(find(im1>.2))=im1(find(im1>.2));
% NWimoverlay(cat(3,im1,im2))

if nargin < 2
    alpha = 0.5;
end

if ndims(Im)~=3
    error('input must be 3d')
end

si = size(Im);
nIm_over = si(3)-1;

if nargin<3 || isempty(ax)
    hfig = figure;
    ax = axes(hfig);
end

imagesc(ax,squeeze(Im(:,:,1))); colormap bone
hold on


% color overlays
if nargin<4 || isempty(colors)
    red = cat(3,ones(si(1:2)),zeros(si(1:2)),zeros(si(1:2)));
    green = cat(3,zeros(si(1:2)),ones(si(1:2)),zeros(si(1:2)));
    yellow = cat(3,ones(si(1:2)),ones(si(1:2)),zeros(si(1:2)));
    orange = cat(3,ones(si(1:2)),0.5*ones(si(1:2)),zeros(si(1:2)));
    blue = cat(3,zeros(si(1:2)),zeros(si(1:2)),ones(si(1:2)));
    fullcol = cat(4,red,green,yellow,orange,blue);
else
    siC = size(colors);
    fullcol = [];
    if siC(1)==3
        for ii=1:siC(2)
            temp = cat(3,colors(1,ii)*ones(si(1:2)),colors(2,ii)*ones(si(1:2)),colors(3,ii)*ones(si(1:2)));
            fullcol = cat(4,fullcol,temp);
        end
    elseif siC(2)==3
        for ii=1:siC(1)
            temp = cat(3,colors(ii,1)*ones(si(1:2)),colors(ii,2)*ones(si(1:2)),colors(ii,3)*ones(si(1:2)));
            fullcol = cat(4,fullcol,temp);
        end
    else
        error('wrong size for colors')
    end
end

h = zeros(1,nIm_over);
for ii=1:nIm_over
    if size(fullcol,4)>1
        h(ii) = imagesc(ax,squeeze(fullcol(:,:,:,mod(ii-1,size(fullcol,4))+1)));
    else
        h(ii) = imagesc(ax,squeeze(fullcol(:,:,:)));
    end
    set(h(ii),'AlphaData',alpha * squeeze(Im(:,:,ii+1)))
end
hold off
