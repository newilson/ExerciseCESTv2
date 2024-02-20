function [p,S] = NWpolyfitim(ndeg,ind_vec,dep_mat)
%
% INPUT:
% ndeg is the polynomial order
% ind_vec is a vector (e.g. 1xt)
% dep_mat is a matrix (e.g. AxBxt)
%
% OUTPUT:
% p are poly coefficients mapped to each dep_mat location (e.g. AxB)
% returned in descending order similar to polyfit.m

npts = length(ind_vec);
si = size(dep_mat);
nd = ndims(dep_mat);

if ~isequal(si(end),npts)
    error('dependent matrix and independent vector do not have equivalent sizes')
end

ind_vec = ind_vec(:); % column vector

opt = 0;
switch opt
    case 0        
        % Vandermonde matrix
        V = zeros(npts,ndeg+1);
        for ii=0:ndeg
            V(:,ndeg+1-ii) = ind_vec.^ii;
        end
        
        % reshape image
        tempdata = reshape(permute(dep_mat,[nd 1:nd-1]),[si(end),numel(dep_mat)/si(end)]);
        
        % solve least squares
        [Q,R] = qr(V,0);
        ptemp = R\(Q'*tempdata); % same as p = V\tempdata;
        
        % reshape coefficients
        p = reshape(ptemp,[ndeg+1,si(1:end-1)]);
        p = permute(p,[2:ndims(p),1]);
        
        if nargout>1
            rtemp = tempdata - V*ptemp;
            r = reshape(rtemp,[npts,si(1:end-1)]);
            S.R = R;
            S.df = max(0,length(ind_vec) - (ndeg+1));
            S.normr = squeeze(sqrt(sum(abs(r).^2,1)));
        end            

    case 1 % do not use; slow, sanity check
        tempdata = reshape(permute(dep_mat,[nd 1:nd-1]),[si(end),numel(dep_mat)/si(end)]);
        p = zeros(prod(si(1:end-1)),ndeg+1); normr = zeros(1,prod(si(1:end-1)));
        for ii=1:prod(si(1:end-1))
            dep_vec = tempdata(:,ii);
            if nargout>1
                [ptemp,Stemp] = polyfit(ind_vec,dep_vec(:),ndeg);
                normr(ii) = Stemp.normr;
            else
                ptemp = polyfit(ind_vec,dep_vec(:),ndeg);
            end
            p(ii,:) = ptemp(:);
        end
        p = reshape(p,[si(1:nd-1),ndeg+1]);
        if nargout>1
            S.R = Stemp.R;
            S.df = Stemp.df;
            S.normr = reshape(normr,si(1:end-1));
        end
end



end
            
