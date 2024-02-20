function [B0maps, hdr2, ampl, phases, slope, intercept, TEs, pathname1, pathname2] = readB0gre(pathname1,pathname2,dofilt)

if nargin<2 || isempty(pathname2)
    [hdr1,ampl,~,pathname1] = readdicomfiles2d;
    [hdr2,phases,~,pathname2] = readdicomfiles2d;
else
    if ~isempty(pathname1)
        [hdr1,ampl] = readdicomfiles2d(pathname1);
    else
        ampl = [];
    end
    [hdr2,phases] = readdicomfiles2d(pathname2);
end
if nargin<3 || isempty(ampl), dofilt = false; end

[nx,ny,nacq] = size(phases);
if ~isempty(ampl)
    ampl = reshape(ampl,[nx,ny,3,nacq/3]);
end
phases = reshape(phases,[nx,ny,3,nacq/3]);
phases = pi * ( (phases-2048.0)/2048.0 );  % Scale phase to be between -pi to +pi
    
if dofilt
    comp = ampl.*exp(1i*phases);
    if length(dofilt)==3
        filt = repmat(NWSiemensRawFilter2d(dofilt(1),dofilt(2),dofilt(3)),[1 1 size(ampl,3)]);
    elseif length(dofilt)==6
        filt = repmat(NWSiemensRawFilter2d(dofilt(1:2),dofilt(3:4),dofilt(5:6)),[1 1 size(ampl,3)]);
    else
        si = size(ampl);
%         filt = repmat(NWSiemensRawFilter2d(size(ampl,1),size(ampl,1),size(ampl,1)/2),[1 1 size(ampl,3)]);
        filt = repmat(NWSiemensRawFilter2d(si(1:2),si(1:2),si(1:2)/2),[1 1 si(3)]);
    end
    comp = fft2c(ifft2c(comp).*filt);
    ampl = abs(comp);
    phases = angle(comp);
end

[B0maps,slope,intercept,TEs] = NWcalcB0gre([],phases,hdr2(1));
B0maps = B0maps/hdr2(1).sf;
% for ii=1:nacq/3
%     B0maps(:,:,ii) = anisodiff(B0maps(:,:,ii),20,50,0.03,1);
% end

return
