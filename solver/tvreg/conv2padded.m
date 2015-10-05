function x = conv2padded(varargin)
%CONV2PADDED  Two-dimensional convolution with padding.
%   Y = CONV2PADDED(X,H) applies 2D filter H to X with constant extension
%   padding. 
%
%   Y = CONV2PADDED(H1,H2,X) first applies 1D filter H1 along the rows and 
%   then applies 1D filter H2 along the columns.
%
%   If X is a 3D array, filtering is done separately on each channel.

% Pascal Getreuer 2009

if nargin == 2       % Function was called as "conv2padded(x,h)"
    x = varargin{1};
    h = varargin{2};
    top = ceil(size(h,1)/2)-1;
    bottom = floor(size(h,1)/2);
    left = ceil(size(h,2)/2)-1;
    right = floor(size(h,2)/2);
elseif nargin == 3   % Function was called as "conv2padded(h1,h2,x)"
    h1 = varargin{1};
    h2 = varargin{2};
    x = varargin{3};
    top = ceil(length(h1)/2)-1;
    bottom = floor(length(h1)/2);
    left = ceil(length(h2)/2)-1;
    right = floor(length(h2)/2);
else
    error('Wrong number of arguments.');
end

% Pad the input image
xPadded = x([ones(1,top),1:size(x,1),size(x,1)+zeros(1,bottom)],...
      [ones(1,left),1:size(x,2),size(x,2)+zeros(1,right)],:);

% Since conv2 cannot handle 3D inputs, we do filtering channel by channel
for p = 1:size(x,3)
    if nargin == 2
        x(:,:,p) = conv2(xPadded(:,:,p),h,'valid');     % Call conv2
    else
        x(:,:,p) = conv2(h1,h2,xPadded(:,:,p),'valid'); % Call conv2
    end
end
