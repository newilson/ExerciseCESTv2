function NWinteractiveCEST(im,xax,im2,xax2)
% 
%
if nargin<3, im2 = []; end
if nargin<4, xax2 = []; end

im = squeeze(im);
si = size(im);
pix = round(si(1:2)/2);
frame = 1;

if ndims(im)>3
    error('input must be 3D')
end

if nargin<2 || isempty(xax)
    xax = 1:si(3);
end

if ~isequal(size(xax),si)
    if ~isequal(length(xax),si(3))
        error('axis length must equal number of frames')
    end
end

if isequal(length(xax),numel(xax))
    xax = repmat(permute(xax(:)',[3 1 2]),[si(1:2) 1]);
end


% Second set (e.g. fits, etc) 
if ~isempty(im2)
    if iscell(im2)
        si2 = size(im2{1});
        nim2 = length(im2);
        for bb=2:nim2
            if size(im2{bb})~=si2
                im2 = [];
                break;
            end
        end
    else
        si2 = size(im2);
    end
    if ~isequal(si(1:2),si2(1:2))
        warning('unmatched sizes')
        im2 = [];
        xax2 = [];
    end
end

if ~isempty(im2) && isempty(xax2)
    xax2 = 1:si2(3);
end

if ~isempty(xax2)
    if isequal(length(xax2),numel(xax2))
        xax2 = repmat(permute(xax2(:)',[3 1 2]),[si2(1:2) 1]);
    end
    if iscell(im2)
        if ~isequal(size(im2{1}),size(xax2))
            warning('mismatched lengths')
            im2 = [];
            xax2 = [];
        end
    else
        if ~isequal(size(im2),size(xax2))
            warning('mismatched lengths')
            im2 = [];
            xax2 = [];
        end
    end
end

% Check if asymmetry curve is possible
if mod(si(3),2)==0 && isequal(-flip(xax(1:si(3)/2)),xax(si(3)/2+1:end))
    allow_asym = true;
    Nim = flip(im(:,:,1:si(3)/2),3);
    Pim = im(:,:,si(3)/2+1:end);
    Aim0 = 100*(Nim-Pim)./repmat(Nim(:,:,end),[1 1 si(3)/2]);
    Aim0 = cat(3,zeros(size(Aim0)),Aim0); % to keep same size as Zspec
    AimN = 100*(Nim-Pim)./Nim;
    AimN = cat(3,zeros(size(AimN)),AimN);
    xlA = [xax(si(3)/2+1) xax(end)];
else
    allow_asym = false;
    Aim0 = [];
    AimN = [];
    xlA = [];
end

% Problems if repeated points in axppm (check this)
for ii=1:si(1)
    for jj=1:si(2)
        reps = find(histc(xax(ii,jj,:),unique(xax(ii,jj,:)))>1);
        xax(ii,jj,reps) = xax(ii,jj,reps)-eps;
        if ~isempty(xax2)
            reps = find(histc(xax2(ii,jj,:),unique(xax2(ii,jj,:)))>1);
            xax2(ii,jj,reps) = xax2(ii,jj,reps)-eps;
        end
    end
end

f = figure('position',[125 90 700 680]);

ax = axes('Parent',f,'position',[.2 .4 .5 .5],'Ydir','reverse');
hIm = imagesc('Parent',ax,'cdata',squeeze(im(:,:,frame))); 

hpoint = impoint(ax,pix);
setColor(hpoint,'r')
fcn = makeConstrainToRectFcn('impoint',[1 si(2)],[1 si(1)]);
setPositionConstraintFcn(hpoint,fcn);
addNewPositionCallback(hpoint,@(pos) pixMove(pos));

cax = caxis(ax);
% caxis(ax,cax);
axis tight, axis equal, axis tight
axis off
colorbar(ax)

ax2 = axes('Parent',f,'position',[.2 .082 .6 .25]);
if isempty(im2)
    hplot = plot(ax2,squeeze(xax(pix(2),pix(1),:)),squeeze(im(pix(2),pix(1),:)),'o');
elseif iscell(im2)
    hold(ax2,'on')
    hplot(1) = plot(ax2,squeeze(xax(pix(2),pix(1),:)),squeeze(im(pix(2),pix(1),:)),'o');
    for jjj=1:nim2
        hplot(jjj+1) = plot(ax2,squeeze(xax2(pix(2),pix(1),:)),squeeze(im2{jjj}(pix(2),pix(1),:)),'-.');
    end
    set(hplot(1),'color',[0.5 0.5 0.5]); % gray
else
    hplot = plot(ax2,squeeze(xax(pix(2),pix(1),:)),squeeze(im(pix(2),pix(1),:)),'o',squeeze(xax2(pix(2),pix(1),:)),squeeze(im2(pix(2),pix(1),:)),'-.');
    set(hplot(1),'color',[0.5 0.5 0.5]); % gray
end
% xlabel('ppm')
axis tight

xl = [min(xax(:))-abs(min(xax(:))/10) max(xax(:))+abs(max(xax(:))/10)];
% yl = [-10 110];
yl = cax;
linelim = [-4096 4096];
hline = imline(ax2,xax(pix(2),pix(1),1)*[1 1],linelim);
setColor(hline,'r')
setPositionConstraintFcn(hline,@(pos) [repmat(mean(pos(:,1)),2,1) linelim(:)]) % vertical
addNewPositionCallback(hline,@(pos) lineMove(pos));


if ndims(im)>2
    s1 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.74 .46 .05 .052],...
        'value',frame,'min',1,'max',si(3),...
        'sliderstep',[1 1]/(si(3)-1),'callback',@nextframe);
    t1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.74 .41 .05 .05],...
        'string',num2str(frame),'fontweight','bold');
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

e3 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.2 .37 .1 .04],...
    'callback',@pixChange);

e4 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .30 .07 .05],...
    'callback',@ylimMin);

e5 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .25 .07 .05],...
    'callback',@ylimMax);

c1 = uicontrol('Parent',f,'style','checkbox','units','normalized','position',[.87 .1 .12 .05],...
    'string','Autoscale','callback',@autosc);

t3 = uicontrol('Parent',f,'style','text','units','normalized','position',[.14 .36 .05 .04],...
    'string','Pixel','HorizontalAlignment','right','fontweight','bold');

te1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .70 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .65 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');

te4 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .30 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .25 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');


bg = uibuttongroup(f,'units','normalized','Position',[.60 .35 .35 .05],'bordertype','none','SelectionChangedFcn',@RadioPlot);
r1 = uicontrol(bg,'Style','radio','String','Set 1','fontweight','bold','units','normalized','position',[.05 .32 .4 .5]);
r2 = uicontrol(bg,'Style','radio','String','Set 2','fontweight','bold','units','normalized','position',[.35 .32 .4 .5]);

if isempty(im2)
   bg.Visible = 'off';
end
   bg.Visible = 'off';


set(f,'visible','on','toolbar','figure')

setPix(pix)
updateStringsIm(cax)
updateStringsPlot(yl)
updateImage
updatePlot

    function nextframe(source,callbackdata)
        frame = round(get(s1,'value'));
        updateImage
        updatePlot
        set(t1,'string',num2str(frame))
        set(s1,'value',frame)
        yl = ylim(ax2);
        temppos = [xax(pix(2),pix(1),frame)*[1;1] linelim(:)];
        setPosition(hline,temppos);
    end

    function caxis_2(source,callbackdata)
        cax = caxis(ax)/2;
        caxis(ax,cax);
        updateStringsIm(cax)
    end
    
    function caxisX2(source,callbackdata)
        cax = caxis(ax)*2;
        caxis(ax,cax);
        updateStringsIm(cax)
    end

    function caxisMax(source,callbackdata)
        cax = caxis(ax);
        if str2double(get(e2,'string'))>cax(1)
            cax(2) = str2double(get(e2,'string'));
            caxis(ax,cax);
        else
            updateStringsIm(cax)
        end
    end

    function caxisMin(source,callbackdata)
        cax = caxis(ax);
        if str2double(get(e1,'string'))<cax(2)
            cax(1) = str2double(get(e1,'string'));
            caxis(ax,cax);
        else
            updateStringsIm(cax)
        end
    end

    function ylimMax(source,callbackdata)
        yl = ylim(ax2);
        if str2double(get(e5,'string'))>yl(1)
            yl(2) = str2double(get(e5,'string'));
            ylim(ax2,yl);
        else
            updateStringsPlot(yl)
        end
    end

    function ylimMin(source,callbackdata)
        yl = ylim(ax2);
        if str2double(get(e4,'string'))<yl(2)
            yl(1) = str2double(get(e4,'string'));
            ylim(ax2,yl);
        else
            updateStringsPlot(yl)
        end
    end

    function updateStringsIm(cax)
        set(e1,'string',num2str(cax(1),'%u'))
        set(e2,'string',num2str(cax(2),'%u'))
    end

    function updateStringsPlot(yl)
        set(e4,'string',num2str(yl(1),'%u'))
        set(e5,'string',num2str(yl(2),'%u'))
        ylim(ax2,yl)
    end

    function updateImage
        if strcmpi(get(get(bg,'SelectedObject'),'String'),'Z spec')
            set(hIm,'cdata',squeeze(im(:,:,get(s1,'Value'))))
        elseif strcmpi(get(get(bg,'SelectedObject'),'String'),'M- Asym')
            set(hIm,'cdata',squeeze(AimN(:,:,abs(get(s1,'Value')-si(3)/2-1)+si(3)/2+1)))
        elseif strcmpi(get(get(bg,'SelectedObject'),'String'),'M0 Asym')
            set(hIm,'cdata',squeeze(Aim0(:,:,abs(get(s1,'Value')-si(3)/2-1)+si(3)/2+1)))
        else
%             warning('unknown radio button')
            set(hIm,'cdata',squeeze(im(:,:,round(get(s1,'Value')))))
        end
    end

    function updatePlot
        set(hplot(1),'xdata',squeeze(xax(pix(2),pix(1),:)),'ydata',squeeze(im(pix(2),pix(1),:)))
        if ~isempty(im2)
            if iscell(im2)
                for iii=1:nim2
                    set(hplot(iii+1),'xdata',squeeze(xax2(pix(2),pix(1),:)),'ydata',squeeze(im2{iii}(pix(2),pix(1),:)))
                end
            else
                set(hplot(2),'xdata',squeeze(xax2(pix(2),pix(1),:)),'ydata',squeeze(im2(pix(2),pix(1),:)))
            end
        end
        xlim(ax2,xl);
%         if strcmpi(get(get(bg,'SelectedObject'),'String'),'Z spec')
%             set(hplot,'xdata',squeeze(axppm(pix(2),pix(1),:)),'ydata',squeeze(im(pix(2),pix(1),:)))
%             xlim(ax2,xl);
%         elseif strcmpi(get(get(bg,'SelectedObject'),'String'),'M- Asym')
%             set(hplot,'ydata',squeeze(AimN(pix(2),pix(1),:)))
%             xlim(ax2,xlA);
%         elseif strcmpi(get(get(bg,'SelectedObject'),'String'),'M0 Asym')
%             set(hplot,'ydata',squeeze(Aim0(pix(2),pix(1),:)))
%             xlim(ax2,xlA);
%         else
%             warning('unknown radio button')
%             set(hplot,'ydata',squeeze(im(pix(2),pix(1),:)))
%             xlim(ax2,xl);
%         end
    end

    function RadioPlot(source,callbackdata)
        updatePlot
        updateImage
    end

    function autosc(source,callbackdata)
        if get(c1,'Value')
            vals = im(pix(2),pix(1),:);
            if ~isempty(im2)
                if iscell(im2)
                    for aa=1:nim2
                        vals = cat(3,vals,im2{aa}(pix(2),pix(1),:));
                    end
                else
                    vals = cat(3,vals,im2(pix(2),pix(1),:));
                end
            end
            vals = vals(:);
            minval = min(vals);
            maxval = max(vals);
            if minval<0
                minyl = 1.1*minval;
            else
                minyl = 0.9*minval;
            end
            if maxval<0
                maxyl = 0.9*maxval;
            else
                maxyl = 1.1*maxval;
            end
            yl = [minyl maxyl];
            updateStringsPlot(yl)
        end
    end
    
    function pixChange(source,callbackdata)
        pixorig = pix;
        value = get(e3,'String');
        expression = ['[' value ']'];
        try
            temp = eval(expression);
            if length(temp)~=2 || any(temp<1) || temp(2)>size(im,1) || temp(1)>size(im,2)
                warndlg('Invalid Pixel')
                setPix(pixorig)
            else
                pix = flip(round(temp));
                setPosition(hpoint,pix)
                updatePlot
                autosc
            end
        catch
            warndlg('Invalid Pixel')
            setPix(pixorig)
        end        
    end

    function setPix(pixel)
        set(e3,'String',[num2str(pixel(2)) ', ' num2str(pixel(1))])
    end

    function pixMove(pos)
        pix = round(pos);
        setPosition(hpoint,pix)
        setPix(pix)
        updatePlot
        autosc
    end

    function lineMove(pos)
        [~,ind] = min(abs(xax(pix(2),pix(1),:) - pos(1)));
        frame = ind;
        set(s1,'Value',frame)
        set(t1,'String',num2str(frame))        
        updateImage
    end    

    function printax(source,callbackdata)
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
    end
  
end
