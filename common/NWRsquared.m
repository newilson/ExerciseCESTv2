function Rsq = NWRsquared(data,fit,dim)
%
% Rsq = NWRsquared(data,fit,dim)
%
% data is the acquired data, fit is the fitted data
% dim is the dimension that was fitted

si = size(data);
if ~isequal(si,size(fit))
    error('acquired and fit data must be the same size')
end
if nargin<3 || isempty(dim)
    dim = 1;
end

SoSres = sum((data-fit).^2,dim);
SoStot = (si(dim)-1)*var(data,0,dim);

Rsq = 1 - SoSres./(SoStot+eps);
end