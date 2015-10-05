% TV denoising demo using tvrestore, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> tvdenoise_demo

clear all;
NoiseLevel = 0.1; 
lambda = 15;

% Simulate a noisy image
uexact = double(imread('einstein.png'))/255;
f = uexact + randn(size(uexact))*NoiseLevel;

% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Denoising');
compareimages(f,'Input',f,'Denoised');
shg;

% Denoise
% initilization lambda
lambda = 0.7079/(NoiseLevel) + 0.6849/(NoiseLevel^2);
for i = 1:10
u = tvdenoise(f,lambda,[],[],[],@tvregsimpleplot);
tmp = f(:)-u(:);
tmp_lambda(i) = lambda;
% lambda update
lambda = lambda*sqrt(sum(tmp.^2)/length(u(:)))/NoiseLevel;
if abs(tmp_lambda(i)-lambda)/tmp_lambda(i) < 1e-2
    break;
end
norm(i) = sum((u(:)-uexact(:)).^2)/(length(u(:)));
end
compareimages(f,'Input',u,'Denoised');
