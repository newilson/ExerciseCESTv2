function NWinteractiveCEST3d(im,xax,im2,xax2,title)
% 
%
if nargin<3, im2 = []; end
if nargin<4, xax2 = []; end
if nargin<5, title = []; end

hsize = 3;
hshape = 0.5;

% im = squeeze(im);
imorig = im;
si = size(im);
pix = round(si(1:2)/2);
frame = 1;
slice = 1;

if ndims(im)>4
    error('input must be 4D')
end

if nargin<2 || isempty(xax)
    xax = 1:si(end);
end

if ~isequal(size(xax),si)
    if ~isequal(length(xax),si(end))
        error('axis length must equal number of frames')
    end
end

if isequal(length(xax),numel(xax))
    xax = repmat(permute(xax(:)',[3 4 1 2]),[si(1:3) 1]);
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
    if ~isequal(si(1:3),si2(1:3))
        warning('unmatched sizes')
        im2 = [];
        xax2 = [];
    end
end

if ~isempty(im2) && isempty(xax2)
    xax2 = 1:si2(end);
end

if ~isempty(xax2)
    if isequal(length(xax2),numel(xax2))
        xax2 = repmat(permute(xax2(:)',[3 4 1 2]),[si2(1:3) 1]);
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

%Check if asymmetry curve is possible
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

%Problems if repeated points in axppm (check this)


for ii=1:si(1)
    for jj=1:si(2)
        for kk=1:si(3)
            reps = find(histc(xax(ii,jj,kk,:),unique(xax(ii,jj,kk,:)))>1);
            xax(ii,jj,kk,reps) = xax(ii,jj,kk,reps)-eps;
            if ~isempty(xax2)
                reps = find(histc(xax2(ii,jj,kk,:),unique(xax2(ii,jj,kk,:)))>1);
                xax2(ii,jj,kk,reps) = xax2(ii,jj,kk,reps)-eps;
            end
        end
    end
end

f = figure('position',[125 90 700 680]);
if ~isempty(title)
    set(f,'Name',title,'NumberTitle','off');
end

ax = axes('Parent',f,'position',[.2 .4 .5 .5],'Ydir','reverse');
hIm = imagesc('Parent',ax,'cdata',squeeze(im(:,:,slice,frame))); 

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
    hplot = plot(ax2,squeeze(xax(pix(2),pix(1),slice,:)),squeeze(im(pix(2),pix(1),slice,:)),'o');
elseif iscell(im2)
    hold(ax2,'on')
    hplot(1) = plot(ax2,squeeze(xax(pix(2),pix(1),slice,:)),squeeze(im(pix(2),pix(1),slice,:)),'o');
    for jjj=1:nim2
        hplot(jjj+1) = plot(ax2,squeeze(xax2(pix(2),pix(1),slice,:)),squeeze(im2{jjj}(pix(2),pix(1),slice,:)),'-.');
    end
    set(hplot(1),'color',[0.5 0.5 0.5]); % gray
else
    hplot = plot(ax2,squeeze(xax(pix(2),pix(1),slice,:)),squeeze(im(pix(2),pix(1),slice,:)),'o',squeeze(xax2(pix(2),pix(1),slice,:)),squeeze(im2(pix(2),pix(1),slice,:)),'-.');
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


if ndims(im)>3 && si(3)>1
    s1 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.74 .46 .05 .052],...
        'value',slice,'min',1,'max',si(3),...
        'sliderstep',[1 1]/(si(3)-1),'callback',@nextslice);
    t1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.74 .41 .05 .05],...
        'string',num2str(slice),'fontweight','bold');
end
b0 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .93 .07 .05],...
    'callback',@printax,'string','Print','fontweight','bold');

b1 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .85 .05 .05],...
    'callback',@caxis_2,'string','/2');

b2 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .80 .05 .05],...
    'callback',@caxisX2,'string','x2');

b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.81 .93 .07 .05],...
    'callback',@exportcsv,'string','Export','fontweight','bold');

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
   

% from NWimoverlayGUI
pop = uicontrol('Parent',f,'style','popup','units','normalized','position',[.85 .50 .08 .05],...
   'String',{'None','Avg','Blur','Sharp'},'Callback', @doFilt);

tpop = uicontrol('Parent',f,'style','text','units','normalized','position',[.85 .55 .05 .05],...
   'string','Filter','HorizontalAlignment','left','fontweight','bold');

es1 = uicontrol('Parent',f,'Style','edit','units','normalized','Position',[.90 .45 .05 .04],...
    'callback',@doFilt,'string',num2str(hsize));

tes1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.85 .45 .05 .04],...
    'string','Size','HorizontalAlignment','left','fontweight','normal');

es2 = uicontrol('Parent',f,'Style','edit','units','normalized','Position',[.90 .40 .05 .04],...
    'callback',@doFilt,'String',num2str(hshape));

tes2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.85 .40 .05 .04],...
    'string','Shape','HorizontalAlignment','left','fontweight','normal');


set(f,'visible','on','toolbar','figure')

setPix(pix)
updateStringsIm(cax)
updateStringsPlot(yl)
updateImage
updatePlot

    function nextslice(source,callbackdata)
        slice = round(get(source,'value'));
        set(t1,'string',num2str(slice))
        set(s1,'value',slice)
        updateImage
        updatePlot
        autosc
    end

    function nextframe(source,callbackdata)
        frame = round(get(s1,'value'));
        updateImage
        updatePlot
        set(t1,'string',num2str(frame))
        set(s1,'value',frame)
        yl = ylim(ax2);
        temppos = [xax(pix(2),pix(1),frame)*[1;1] linelim(:)];
        setPosition(hline,temppos);
        autosc
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

    function doFilt(source,callbackdata) % from NWimoverlayGUI
        hsize = str2double(es1.String);
        hshape = str2double(es2.String);
        if isequal(pop.Value,2)
            filt = fspecial('average',hsize);
            im = imfilter(imorig,filt,'replicate');
        elseif isequal(pop.Value,3)
            im = imgaussfilt(imorig,hshape,'FilterSize',hsize,'Padding','replicate');     
        elseif isequal(pop.Value,4)
            for ii=1:size(im,4)
                for jj=1:size(im,3)
                    im(:,:,jj,ii) = imsharpen(imorig(:,:,jj,ii),'Radius',hsize/2,'Amount',hshape);
                end
            end
        elseif isequal(pop.Value,1)
            im = imorig;
        end
        
        slice = get(s1,'Value');
        updateImage
        updatePlot
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
            set(hIm,'cdata',squeeze(im(:,:,slice,frame)))
        end
    end

    function updatePlot
        set(hplot(1),'xdata',squeeze(xax(pix(2),pix(1),slice,:)),'ydata',squeeze(im(pix(2),pix(1),slice,:)))
        if ~isempty(im2)
            if iscell(im2)
                for iii=1:nim2
                    set(hplot(iii+1),'xdata',squeeze(xax2(pix(2),pix(1),slice,:)),'ydata',squeeze(im2{iii}(pix(2),pix(1),slice,:)))
                end
            else
                set(hplot(2),'xdata',squeeze(xax2(pix(2),pix(1),slice,:)),'ydata',squeeze(im2(pix(2),pix(1),slice,:)))
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
            vals = im(pix(2),pix(1),slice,:);
            if ~isempty(im2)
                if iscell(im2)
                    for aa=1:nim2
                        vals = cat(4,vals,im2{aa}(pix(2),pix(1),slice,:));
                    end
                else
                    vals = cat(4,vals,im2(pix(2),pix(1),slice,:));
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
        [~,ind] = min(abs(xax(pix(2),pix(1),slice,:) - pos(1)));
        frame = ind;
%         set(s1,'Value',frame)
%         set(t1,'String',num2str(frame))        
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

    function exportcsv(source,callbackdata)
        if ~isempty(title)
            fname = [title '_pix' num2str(pix(2)) '_' num2str(pix(1)) '_' num2str(slice) '.csv'];
        else
            fname = ['pix' num2str(pix(2)) '_' num2str(pix(1)) '_' num2str(slice) '.csv'];
        end
        col1 = squeeze(xax(pix(2),pix(1),slice,:)); col1 = col1(:);
        col2 = squeeze(im(pix(2),pix(1),slice,:)); col2 = col2(:);
        
        if ~isempty(im2)
            col3 = squeeze(xax2(pix(2),pix(1),slice,:)); col3 = col3(:);
            col4 = im2(pix(2),pix(1),slice,:); col4 = col4(:);
            VariableNames = {'xax','sig','xax2','sig2'};
            T = table(col1,col2,col3,col4,'VariableNames',VariableNames);
        else
            VariableNames = {'xax','sig'};
            T = table(col1,col2,'VariableNames',VariableNames);
        end
        writetable(T,fname,'Delimite',',');
    end
  
end
