function [slope,intercept,Rsq,fit] = NWlinfit_series(ind_vec,dep_mat)
%
% [slope,intercept,Rsq,fit] = NWlinfit_series(ind_vec,dep_mat)
%
% INPUT:
% ind_vec is a vector (e.g. 1xt)
% dep_mat is a matrix (e.g. AxBxt)
%
% OUTPUT:
% slope and intercept mapped to each dep_mat location (e.g. AxB)
% Rsq is the coefficient of determination

n = length(ind_vec);
si = size(dep_mat);
nd = ndims(dep_mat);

if ~isequal(si(end),n)
    error('last dimension of dep_mat must have same length as ind_vec')
end

ind_vec = ind_vec(:); % column vector

opt = 0;
switch opt
    case 0 % fast way
        A = [ones(n,1),ind_vec];
        temp = reshape(permute(dep_mat,[nd 1:nd-1]),[si(end),numel(dep_mat)/si(end)]);
        coef = A \ temp;
        slope = reshape(coef(2,:),si(1:nd-1));
        intercept = reshape(coef(1,:),si(1:nd-1));
        
    case 1 % slow way
        temp = reshape(permute(dep_mat,[nd 1:nd-1]),[si(end),numel(dep_mat)/si(end)]);
        slope = zeros(1,prod(si(1:end-1))); intercept = 0*slope;
        for ii=1:prod(si(1:end-1))
            dep = temp(:,ii);
            p = polyfit(ind_vec,dep(:),1);
            slope(ii) = p(1);
            intercept(ii) = p(2);
        end
        slope = reshape(slope,si(1:nd-1));
        intercept = reshape(intercept,si(1:nd-1));
        
    case 2 % very slow way, only works for 3d images right now
        if nd==4
            slope = zeros(si(1:nd-1));
            intercept = zeros(si(1:nd-1));
            for ii=1:si(1)
                for jj=1:si(2)
                    for kk=1:si(3)
                        dep = dep_mat(ii,jj,kk,:);
                        p = polyfit(ind_vec,dep(:),1);
                        slope(ii,jj,kk) = p(1);
                        intercept(ii,jj,kk) = p(2);
                    end
                end
            end
        else
            warning('not yet implemented')
        end
end

if nargout>2
    fit = zeros(size(dep_mat));
    indices = repmat({':'},[1,nd-1]); % unknown number of colon operators
    for ii=1:n
        tempfit = slope*ind_vec(ii) + intercept;
        fit(indices{:},ii) = tempfit;
    end
    Rsq = NWRsquared(dep_mat,fit,nd);
end
    
end