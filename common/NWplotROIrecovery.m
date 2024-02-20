function colors = NWplotROIrecovery(timepoints,fitdata,means,stds,ref,rois,colors,appstruct)


h = figure('Name','ROIs ','NumberTitle','off');
ax = axes(h); hold on

if nargin<8, mainDir = pwd; else, mainDir = appstruct.mainDir; end
if nargin<6, rois = []; end
if nargin<5, ref = []; end
if nargin<4, stds = 0*means; end

if isempty(fitdata)
    nrois = size(means,1);
else
    nrois = size(fitdata,2);
end

if nargin<7 || isempty(colors) || ~isequal(size(colors,1),nrois), colors = get(gca,'ColorOrder'); end

plotmax = 1.2*max(means(:)+stds(:));

for ii=1:nrois
    if ~isempty(fitdata)
        plot(ax,timepoints,fitdata(:,ii),'color',colors(ii,:),'linewidth',1.5)
    end
%     errorbar(ax,timepoints,means(ii,:),stds(ii,:),'marker','o','color',colors(ii,:))
    plot(ax,timepoints,means(ii,:),'marker','o','color',colors(ii,:))
    ylim([0 plotmax])
    xlim([timepoints(1)-10 timepoints(end)+10])
    grid on
    xlabel('Time (s)')
    ylabel('% Asym')
end
ax2 = axes(h,'Position',[0.67 0.65 0.25 0.25]);
NWimoverlay(cat(3,ref,rois),0.5,ax2,colors);
colormap('bone')
axis off
axis equal

b0 = uicontrol('Parent',h,'style','pushbutton','units','normalized','position',[.9 .93 .07 .05],...
    'callback',{@printax,mainDir},'string','Print','fontweight','bold');

e1 = uicontrol('Parent',h,'style','edit','units','normalized','position',[.21 .93 .05 .05],...
    'callback',@ylimMin);

e2 = uicontrol('Parent',h,'style','edit','units','normalized','position',[.34 .93 .05 .05],...
    'callback',@ylimMax);

te1 = uicontrol('Parent',h,'style','text','units','normalized','position',[.15 .93 .06 .05],...
    'string','YMin','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',h,'style','text','units','normalized','position',[.28 .93 .06 .05],...
    'string','YMax','HorizontalAlignment','right','fontweight','bold');

updateStrings(ylim(ax))


function printax(source,callbackdata,mainDir)
    prompt = {'FullName (no extension)','Format (vector: eps,pdf / bitmap: tiff,png,bmp,jpeg)','DPI','Renderer (painters or opengl)'};
    defvals = {'recovery','png','300','painters'};
    nlines = 1;
    vals = inputdlg(prompt,'Options',nlines,defvals);
    if ~isempty(vals)
        fname = vals{1};
        if strcmp(vals{2},'eps')
            form = '-depsc';
        else
            form = ['-d' vals{2}];
        end
        res = ['-r' vals{3}];
        rend = ['-' vals{4}];
        print(h,fullfile(mainDir,fname),form,res,rend,'-noui')
    end
end

function ylimMax(source,callbackdata)
    yl = ylim(ax);
    if str2double(get(e2,'string'))>yl(1)
        yl(2) = str2double(get(e2,'string'));
        ylim(ax,yl);
    else
        updateStrings(yl)
    end
end

function ylimMin(source,callbackdata)
    yl = ylim(ax);
    if str2double(get(e1,'string'))<yl(2)
        yl(1) = str2double(get(e1,'string'));
        ylim(ax,yl);
    else
        updateStrings(yl)
    end
end

function updateStrings(yl)
    set(e1,'string',num2str(yl(1),3))
    set(e2,'string',num2str(yl(2),3))
end

end