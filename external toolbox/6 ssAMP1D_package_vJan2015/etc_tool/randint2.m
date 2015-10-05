function out = randint2(varargin)
%
%
%WARNING: This is an obsolete function and may be removed in the future.
%         Please use RANDI instead.
%
%
%RANDINT Generate matrix of uniformly distributed random integers.
%   OUT = RANDINT generates a "0" or "1" with equal probability.
%
%   OUT = RANDINT(M) generates an M-by-M matrix of random binary numbers.
%   "0" and "1" occur with equal probability.
%
%   OUT = RANDINT(M,N) generates an M-by-N matrix of random binary numbers.
%   "0" and "1" occur with equal probability.
%
%   OUT = RANDINT(M,N,IRANGE) generates an M-by-N matrix of random integers.
%
%   IRANGE can be either a scalar or a two-element vector:
%   Scalar : If IRANGE is a positive integer, then the output integer
%            range is [0, IRANGE-1].  If IRANGE is a negative integer,
%            then the output integer range is [IRANGE+1, 0].
%   Vector : If IRANGE is a two-element vector, then the output
%            integer range is [IRANGE(1), IRANGE(2)].
%
%   OUT = RANDINT(M,N,IRANGE,STATE) causes RAND to use the generator
%   determined by the 'state' method, and initializes the state of that
%   generator using the value of STATE.
%
%   Examples:
%       r1 = randint(2,3)                 
%       r2 = randint(2,3,4)
%       r3 = randint(2,3,-4)              
%       r4 = randint(2,3,[-2 2])
%
%   See also RAND, RANDSRC, RANDERR.

%   Copyright 1996-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/08/29 23:30:29 $


% Basic function setup.
error(nargchk(0,4,nargin,'struct'));

% --- Placeholder for the signature string.
sigStr = '';
m = [];
n = [];
range = [];
state = [];

% --- Identify string and numeric arguments
for i=1:nargin
   if(i>1)
      sigStr(size(sigStr,2)+1) = '/';
   end;
   % --- Assign the string and numeric flags
   if(isnumeric(varargin{i}))
      sigStr(size(sigStr,2)+1) = 'n';
   else
      error('comm:randint:InvalidArg','Only numeric arguments are accepted.');
   end;
end;

% --- Identify parameter signatures and assign values to variables
switch sigStr
   % --- randint
   case ''

   % --- randint(m)
   case 'n'
      m		= varargin{1};

	% --- randint(m, n)
	case 'n/n'
      m		= varargin{1};
      n		= varargin{2};

	% --- randint(m, n, range)
	case 'n/n/n'
      m		= varargin{1};
      n  	= varargin{2};
      range = varargin{3};

	% --- randint(m, n, range, state)
	case 'n/n/n/n'
      m		= varargin{1};
      n		= varargin{2};
      range = varargin{3};
      state = varargin{4};

   % --- If the parameter list does not match one of these signatures.
   otherwise
      error('comm:randint:InvalidSyntax','Syntax error.');
end;

if isempty(m)
   m = 1;
end
if isempty(n)
   n = m;
end
if isempty(range)
   range = [0, 1];
end

len_range = size(range,1) * size(range,2);

% Typical error-checking.
if all(length(m) > 1) || all(length(n) > 1)
   error('comm:randint:InvalidMatrixDims','Matrix dimensions must be scalars.');
elseif (floor(m) ~= m) || (floor(n) ~= n) || (~isreal(m)) || (~isreal(n))
   error('comm:randint:NonIntegerMatrixDims','Matrix dimensions must be real integers.');
elseif (m < 0) || (n < 0)
   error('comm:randint:NonPositiveMatrixDims','Matrix dimensions must be positive.');
elseif (~isfinite(m)) || (~isfinite(n))
   error('comm:randint:NonFiniteMatrixDims','Matrix dimensions must be finite.');
elseif len_range > 2
   error('comm:randint:InvalidIrange','The IRANGE parameter should contain no more than two elements.');
elseif max(max(floor(range) ~= range)) || (~isreal(range)) || all(~isfinite(range))
   error('comm:randint:NonIntIrange','The IRANGE parameter must only contain real finite integers.');
end

% If the IRANGE is specified as a scalar.
if len_range < 2
    if range < 0
        range = [range+1, 0];
    elseif range > 0
        range = [0, range-1];
    else
        range = [0, 0];    % Special case of zero range.
    end
end

% Make sure IRANGE is ordered properly.
range = sort(range);

% Calculate the range the distance for the random number generator.
distance = range(2) - range(1);

% Set the initial state if specified.
if ~isempty(state)
   rand('state', state);
end

% Generate the random numbers.
r = floor(rand(m, n) * (distance+1));

% Offset the numbers to the specified value.
out = ones(m,n)*range(1);
out = out + r;

% [EOF] randint.m
