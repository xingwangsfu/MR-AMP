function compareimages(A,ATitle,B,BTitle)
%COMPAREIMAGES   Displays two images side by side with linked axes
%   COMPAREIMAGES(A,B) displays images A and B, where A and B are either
%   grayscale or RGB color images with values in [0,1].  The images are
%   displayed with linked axes for convenient panning and zooming.
%
%   COMPAREIMAGES(A,'A title',B,'B title') specifies titles above the
%   images.
%
%   See also linkaxes.

% Pascal Getreuer 2009

if nargin == 2
    B = ATitle;
    ATitle = '';
    BTitle = '';
elseif nargin ~= 4
    error('Must have 2 or 4 input arguments.');
end

ax(1) = subplot(1,2,1);
hold off

if ndims(A) == 2
    image(A*255);
    colormap(gray(256));
elseif ndims(A) == 3
    image(min(max(A,0),1));
end

set(gca,'Units','Normalized','Position',[0,0.1,0.5,0.8]);
axis image
axis off
title(ATitle);
zoom;

ax(2) = subplot(1,2,2);
hold off

if ndims(B) == 2
    image(B*255);
    colormap(gray(256));
elseif ndims(B) == 3
    image(min(max(B,0),1));
end

set(gca,'Units','Normalized','Position',[0.5,0.1,0.5,0.8]);
axis image
axis off
title(BTitle);

linkaxes(ax,'xy');
