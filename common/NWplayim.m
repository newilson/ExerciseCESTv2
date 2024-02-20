function f = NWplayim(im,cmap,showcbar,fig_title,aspectratio)

im = squeeze(im);
si = size(im);

if ndims(im)>3
    im = reshape(im,[si(1:2), numel(im)/prod(si(1:2))]);
%     error('input must be 3D')
end
si = size(im); % update

% f = figure('position',[1921 -662 1080 1834]);
% f = figure('position',[2561 -789 1080 1782]);
f = figure('position',[200 200 800 450]);
if nargin>=4 && ~isempty(fig_title)
    set(f,'Name',fig_title,'NumberTitle','off');
end
if nargin<5, aspectratio = []; end
% ax = axes('Parent',f,'position',[.1 .2 .8 .8*si(1)/si(2)]);
ax = axes('Parent',f,'position',[.1 .1 .7 .8]);
h = imagesc('Parent',ax,'cdata',squeeze(im(:,:,1))); 
cax = [1.1*min(im(:)), 0.9*max(im(:))];
cax = sort(cax); % in case min>max
caxis(ax,cax);
set(ax,'Ydir','reverse')
if ~isempty(aspectratio)
    daspect(aspectratio)
else
    axis tight, axis equal, axis tight
end
axis off
if nargin<2 || isempty(cmap)
    colormap('bone')
else
    try
        colormap(cmap)
    catch
        cmap = 'bone';
        colormap(cmap)
    end
end
if nargin>2 && ~isempty(showcbar) && showcbar
    colorbar(ax)
end

if ndims(im)>2
    s1 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.9 .1 .05 .55],...
        'value',1,'min',1,'max',si(3),...
        'sliderstep',[1 1]/(si(3)-1),'callback',@nextslice);
    t1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.9 .05 .05 .05],...
        'string',num2str(1));
end
b0 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .93 .07 .05],...
    'callback',@printax,'string','Print','fontweight','bold');

b1 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .85 .05 .05],...
    'callback',@caxis_2,'string','/2');

b2 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .80 .05 .05],...
    'callback',@caxisX2,'string','x2');

e1 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .70 .07 .05],...
    'callback',@caxisMin);

e2 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .65 .07 .05],...
    'callback',@caxisMax);
updateStrings(cax)

te1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .70 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .65 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');

set(f,'visible','on','toolbar','figure')

    function nextslice(source,callbackdata)
        slice = round(get(source,'value'));
        set(h,'cdata',squeeze(im(:,:,slice)))
        set(t1,'string',num2str(slice))
        set(s1,'value',slice)
    end

    function caxis_2(source,callbackdata)
        cax = caxis(ax)/2;
        caxis(ax,cax);
        updateStrings(cax)
    end
    
    function caxisX2(source,callbackdata)
        cax = caxis(ax)*2;
        caxis(ax,cax);
        updateStrings(cax)
    end

    function caxisMax(source,callbackdata)
        cax = caxis(ax);
        if str2double(get(e2,'string'))>cax(1)
            cax(2) = str2double(get(e2,'string'));
            caxis(ax,cax);
        else
            updateStrings(cax)
        end
    end

    function caxisMin(source,callbackdata)
        cax = caxis(ax);
        if str2double(get(e1,'string'))<cax(2)
            cax(1) = str2double(get(e1,'string'));
            caxis(ax,cax);
        else
            updateStrings(cax)
        end
    end

    function updateStrings(cax)
        set(e1,'string',num2str(cax(1),3))
        set(e2,'string',num2str(cax(2),3))
    end

    function printax(source,callbackdata)
%         cax = caxis(ax);
%         f2 = figure('visible','off');
%         ax2 = axes('Parent',f2);
%         copyobj(allchild(ax),ax2);
%         
%         set(ax2,'Ydir','reverse')
%         set(ax2,'LooseInset',get(ax2,'TightInset'))
%         axis(ax2,'tight'), axis(ax2,'equal'), axis(ax2,'tight')
%         axis(ax2,'off')
%         colormap(ax2,cmap)
%         caxis(ax2,cax);
        prompt = {'FullName (no extension)','Format (vector: eps,pdf / bitmap: tiff,png,bmp,jpeg)','DPI','Renderer (painters or opengl)'};
        defvals = {fullfile(pwd,'myFigure'),'eps','300','painters'};
        nlines = 1;
        vals = inputdlg(prompt,'Options',nlines,defvals);
        if ~isempty(vals)
            fname = vals{1};
            if strcmp(vals{2},'eps') && ~strcmp(cmap,'bone')
                form = '-depsc';
            else
                form = ['-d' vals{2}];
            end
            res = ['-r' vals{3}];
            rend = ['-' vals{4}]; 
            print(f,fname,form,res,rend,'-noui')
        end
%         close(f2);
    end
  
end