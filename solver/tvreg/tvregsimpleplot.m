function tvregsimpleplot(state, iter, delta, u, h) %#ok<INUSL>
%TVREGSIMPLEPLOT  A simple plot callback to use with tvreg
%   TVREGSIMPLEPLOT(state, iter, delta, u) plots the grayscale or color 
%   image u on the current axis, and displays the iteration number iter in
%   the axis title.  This function can be passed as a plot callback to 
%   tvdenoise, tvdeconv, tvinpaint, tvrestore, or tvreg.
%
%   TVREGSIMPLEPLOT(...,h) specifies the axis handle to use.  To pass the
%   handle to tvdenoise, etc., do for example
%
%       h = gca;  % Save the axis handle
%       plotfun = @(state,iter,delta,u) ...
%                     tvregsimpleplot(state,iter,delta,u,h);
%       tvdenoise(u,lambda,[],[],[],plotfun);
%
%   See also tvreg.

% Pascal Getreuer 2010

if nargin < 5
    h = gca;
end

if size(u,3) == 3
    % Display a color image
    image(min(max(u,0),1), 'Parent', h);
else
    % Display a grayscale image
    image(u*255, 'Parent', h);
    colormap(gray(256));
end

axis(h, 'image');
axis(h, 'off');
title(h, sprintf('iter %d', iter));
