% TV inpainting demo using tvrestore, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> tvinpaint_demo


lambda = 1e4;

uexact = double(imread('einstein.png'))/255;

% Construct inpainting region
[x,y] = meshgrid(1:size(uexact,2),1:size(uexact,1));
[th,r] = cart2pol(x-size(uexact,2)/2,y-size(uexact,1)/2);
D = (sin(r/2+th) > 0.75);

f = uexact;
f(D) = rand(nnz(D),1);

% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Inpainting');
compareimages(f,'Input',f,'Inpainted');
shg;

% Inpaint
u = tvinpaint(f,lambda,D,[],[],[],@tvregsimpleplot);

title('Inpainted');
