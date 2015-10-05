function chanvesesimpleplot(state, iter, delta, phi, h) %#ok<INUSL>
%CHANVESESIMPLEPLOT  A simple plot callback to use with chanvese
%   CHANVESESIMPLEPLOT(state, iter, delta, phi) plots the segmentation
%   described by phi on the current axis, and displays the current
%   iteration number iter in the axis title.  This function can be used as
%   a plot callback with chanvese.
%
%   CHANVESESIMPLEPLOT(...,h) specifies the axis handle to use.  To pass 
%   the handle to chanvese, do for example
%
%       h = gca;  % Save the axis handle
%       opt.plotfun = @(state,iter,delta,phi) ...
%                     chanvesesimpleplot(state,iter,delta,phi,h);
%       chanvese(f,[],opt);
%
%   See also chanvese.

% Pascal Getreuer 2010


if nargin < 5
    h = gca;
end

imagesc(phi >= 0, 'Parent', h);
axis(h, 'image');
axis(h, 'off');
title(h, sprintf('iter %d', iter));
