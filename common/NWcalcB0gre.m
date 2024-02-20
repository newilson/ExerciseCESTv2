function [map,slope,intercept,TE,Rsq,fitted,phases] = NWcalcB0gre(mag,phases,hdr,mask_flag,opt)

slope = 0; intercept = 0;
if nargin<4, mask_flag = false; end
if nargin<5, opt = 3; end

[nx,ny,ne,nt] = size(phases);
% if ne~=3, error('incorrect sizes of phases'), end
if ~isempty(mag) && size(mag)~=size(phases), error('mag and phases must be same size'), end

if length(hdr)>1
    TE = zeros(1,length(hdr));
    for ii=1:length(hdr)
        TE(ii) = hdr{ii}.TEms/1000;
    end
else
    TE = hdr.TEms/1000;
end
TE = TE(1:ne);
delta_te = diff(TE);

mask = ones(nx,ny,nt);
map12 = zeros(nx,ny,nt); map23 = zeros(nx,ny,nt); map13 = zeros(nx,ny,nt);
map1 = zeros(nx,ny,nt); map2 = zeros(nx,ny,nt); map3 = zeros(nx,ny,nt);
mapf = zeros(nx,ny,nt);
if opt<3
    for ii=1:nt
        if mask_flag
            mask(:,:,ii) = squeeze(mag(:,:,1,ii))>mean(colNW(mag(:,:,1,ii)));
        end
        switch opt
            case 0
                map12(:,:,ii) = mask(:,:,ii).*mexunwrap2DSRNP(wrapToPi(squeeze(phases(:,:,2,ii)-phases(:,:,1,ii))))/(delta_te(1)*2*pi);
                map23(:,:,ii) = mask(:,:,ii).*mexunwrap2DSRNP(wrapToPi(squeeze(phases(:,:,3,ii)-phases(:,:,2,ii))))/(delta_te(2)*2*pi);
                map13(:,:,ii) = mask(:,:,ii).*mexunwrap2DSRNP(wrapToPi(squeeze(phases(:,:,3,ii)-phases(:,:,1,ii))))/(sum(delta_te(:))*2*pi);
            case 1
                map12(:,:,ii) = mask(:,:,ii).*get2DUnwrap(wrapToPi(squeeze(phases(:,:,2,ii)-phases(:,:,1,ii))),squeeze(mag(:,:,1,ii)))/(delta_te(1)*2*pi);
                map23(:,:,ii) = mask(:,:,ii).*get2DUnwrap(wrapToPi(squeeze(phases(:,:,3,ii)-phases(:,:,2,ii))),squeeze(mag(:,:,2,ii)))/(delta_te(2)*2*pi);
                map13(:,:,ii) = mask(:,:,ii).*get2DUnwrap(wrapToPi(squeeze(phases(:,:,3,ii)-phases(:,:,1,ii))),squeeze(mag(:,:,3,ii)))/(sum(delta_te(:))*2*pi);
            case 2
                for jj=1:nx
                    for kk=1:ny
                        p = polyfit(TE(:),colNW(phases(jj,kk,:,ii)),1);
                        mapf(jj,kk,ii) = p(1)/(2*pi);
                    end
                end
        end
    end
else
    phases = unwrap(phases,pi,3);
    if mask_flag
        mask = squeeze(mag(:,:,1,:)>mean(colNW(mag(:,:,1,:))));
    end
    if nargout<4
        [slope,intercept] = NWlinfit_series(TE(:),permute(phases,[1 2 4 3]));
    else
        [slope,intercept,Rsq,fitted] = NWlinfit_series(TE(:),permute(phases,[1 2 4 3]));
    end
    mapf = slope/(2*pi);
end

if opt<2
    map = cat(4,map12,map23,map13);
    % map = squeeze(mean(map,4));
else
    map = mask.*mapf;
end