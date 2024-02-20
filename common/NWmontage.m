function montIm = NWmontage(Im,mSize)

% based on Matlab function montage.m
%
%   Input
%   ----------
%   Im             >3D matrix only (unlike montage.m)
%
%   mSize          A 2-element vector, [NROWS NCOLS], specifying the number
%                  of rows and columns in the montage. Use NaNs to have 
%                  MONTAGE calculate the size in a particular dimension in
%                  a way that includes all the images in the montage. For
%                  example, if 'Size' is [2 NaN], MONTAGE creates a montage
%                  with 2 rows and the number of columns necessary to
%                  include all of the images.  MONTAGE displays the images
%                  horizontally across columns.
%
%                  Default: MONTAGE calculates the rows and columns so the
%                  images in the montage roughly form a square.

si = size(Im);

if ndims(Im)<3
    error('Image must have >3 dimensions')
elseif ndims(Im)>3
    Im = reshape(Im,si(1),si(2),[]);
end

nRows = size(Im,1);
nCols = size(Im,2);
nFrames = size(Im,3);

if nargin<2 || isempty(mSize) || all(isnan(mSize))
    %Calculate montageSize for the user
    
    % Estimate nMontageColumns and nMontageRows given the desired
    % ratio of Columns to Rows to be one (square montage).
    aspectRatio = 1;
    montageCols = sqrt(aspectRatio * nRows * nFrames / nCols);
    
    % Make sure montage rows and columns are integers. The order in
    % the adjustment matters because the montage image is created
    % horizontally across columns.
    montageCols = ceil(montageCols);
    montageRows = ceil(nFrames / montageCols);
    montageSize = [montageRows montageCols];
    
elseif any(isnan(mSize))
    montageSize = mSize;
    nanIdx = isnan(mSize);
    montageSize(nanIdx) = ceil(nFrames / mSize(~nanIdx));
    
elseif prod(mSize) < nFrames
    error(message('images:montage:sizeTooSmall'));
    
else
    montageSize = mSize;
end


nMontageRows = montageSize(1);
nMontageCols = montageSize(2);

sizeOfBigImage = [nMontageRows*nRows nMontageCols*nCols];
if islogical(Im)
    montIm = false(sizeOfBigImage);
else
    montIm = zeros(sizeOfBigImage,'like',Im);
end

rows = 1 : nRows;
cols = 1 : nCols;
k = 1;

for i = 0 : nMontageRows-1
    for j = 0 : nMontageCols-1,
        if k <= nFrames
            montIm(rows + i * nRows, cols + j * nCols, :) = ...
                Im(:,:,k);
        else
            return;
        end
        k = k + 1;
    end
end
