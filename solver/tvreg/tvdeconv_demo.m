% TV deconvolution demo using tvrestore, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> tvdeconv_demo


BlurRadius = 3;
NoiseLevel = 0.005; 
lambda = 6e3;

uexact = double(imread('einstein.png'))/255;

% Simulate a noisy and blurry image
% Construct a disk shaped blur filter
[x,y] = meshgrid(-BlurRadius:BlurRadius,-BlurRadius:BlurRadius);
K = double(x.^2 + y.^2 <= BlurRadius^2);
K = K/sum(K(:));

% Blur the image
f = conv2padded(uexact,K);

% Add Gaussian noise
f = f + randn(size(f))*NoiseLevel;

% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Deconvolution');
compareimages(f,'Input',f,'Deconvolved');
shg;

% Deblur
u = tvdeconv(f,lambda,K,[],[],[],@tvregsimpleplot);

title('Deconvolved');


