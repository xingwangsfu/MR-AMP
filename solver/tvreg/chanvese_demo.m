% Chan-Vese two-phase segmentation demo, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> chanvese_demo


clear opt;

% Length penalty parameter
opt.mu = 0.18;
% Area penalty
opt.nu = 0;
% Inside fit penalty 
opt.lambda1 = 1;
% Outside fit penalty
opt.lambda2 = 1;

% Convergence tolerance
opt.tol = 1e-4;
% Maximum number of iterations
opt.maxiter = 500;
% Timestep parameter
opt.dt = 0.5;


f = double(imread('wrench.bmp'))/255;

% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','Chan-Vese Segmentation');
compareimages(f,'Input',f,'Segmentation');
shg;

% Variable h is to remember which axis we are plotting on.
% In case the user clicks on other figures during the chanvese computation,
% this makes sure the plotting will still be on this figure.
h = gca;
opt.plotfun = @(state,iter,delta,phi) ...
    chanvesesimpleplot(state,iter,delta,phi,h);

phi = chanvese(f,[],opt);

title(h,'Segmentation');
