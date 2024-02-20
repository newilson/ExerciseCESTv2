function [Y,delta] = NWpolyvalim(p,X,S)
%
% p is [X x Y x Z x ndeg] or [1 x ndeg]
% X is [X x Y x Z x npts] or [1 x npts]
% Y is [X x Y x Z x npts]

npts = 1;

if isvector(p)
    p = p(:);
    si = size(X);
    if si(end)>1
        npts = si(end);
        imsize = si(1:end-1);
    else
        imsize = si;
    end
    p = repmat(p,[1 si]);
    p = permute(p,[2:ndims(p) 1]);
end

if isvector(X)
    X = X(:);
    npts = length(X);
    si = size(p);
    ndeg = si(end);
    imsize = si(1:end-1);
    X = repmat(X,[1 si(1:end-1)]);
    X = permute(X,[2:ndims(X) 1]);
end

if npts>1
    p = repmat(p,[ones(1,ndims(p)) npts]);
    p = permute(p,[1:ndims(p)-2 ndims(p) ndims(p)-1]);
end

si = size(p);
ndeg = si(end)-1;

p = permute(p,[ndims(p) 1:ndims(p)-1]);

if ~isequal(size(X),si(1:end-1))
    error('unequal matrix sizes')
end

% Horner's method
Y = 0*X;
for ii=1:ndeg+1
    Y = Y.*X + squeeze(p(ii,:,:,:,:));
end

if nargout>1
    if nargin<3 || isempty(S)
        warning('S required')
        return
    end
    
    % Vandermonde matrix
    X = X(:);
    V = zeros(length(X),ndeg+1);
    for ii=0:ndeg
        V(:,ndeg+1-ii) = X.^ii;
    end
    
    E = V/S.R;
    e = sqrt(1+sum(E.*E,2));
    
    if S.df>0
        delta = S.normr/sqrt(S.df)*e;
    else
        delta = [];
    end
end

end