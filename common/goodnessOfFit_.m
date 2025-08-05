function fit = goodnessOfFit_(x,xref,measure)

% reimplementation of Matlab function goodnessOfFit
% Seems to be an error by the number of samples for MSE in Matlab version
% NW 01/2025

if ~isvector(x)
    warning('use vectors only')
    x = x(:);
    xref = xref(:);
end

if size(x)~=size(xref)
    error('unequal sizes')
end

if strcmpi(measure,'mse')
    fit = norm(x-xref)^2/length(x);
else
    n1 = norm(x-xref);
    n2 = norm(xref - mean(xref(:)));
    if n2==0 && n1==0
        rat = 0;
    else
        rat = n1/n2;
    end
    if strcmpi(measure,'nmse')
        rat = rat^2;
    end
    % fit = 1- rat; % old way pre 9.8
    fit = rat;
end
