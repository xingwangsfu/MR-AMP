function out = randerr(varargin)
%RANDERR Generate bit error patterns.
%   OUT = RANDERR(M) generates an M-by-M binary matrix with exactly
%   one "1" randomly placed in each row.
%
%   OUT = RANDERR(M,N) generates an M-by-N binary matrix with exactly
%   one "1" randomly placed in each row.
%
%   OUT = RANDERR(M,N,ERRORS) generates an M-by-N binary matrix, with
%   the number of "1"s in any given row determined by ERRORS.
%
%   ERRORS can be either a scalar, row vector or two-row matrix:
%   Scalar  : If ERRORS is a scalar the number of "1"s to place in each
%             row is defined by ERRORS.
%   Vector  : If ERRORS is a row vector, then its elements specify how
%             many "1"s are possible.  Every number of "1"s included in
%             this vector occur with equal probability.
%   Two-Row : If ERRORS is a two-row matrix, the first row specifies how
%   Matrix    many "1"s are possible and the second row specifies the
%             probabilities of each corresponding number of "1"s.  The
%             elements in the second row of ERRORS must sum to one.
%
%   OUT = RANDERR(M,N,ERRORS,S) causes RAND to use the random stream S.  See
%   RandStream for more details.
%
%   OUT = RANDERR(M,N,ERRORS,STATE) causes RAND to use the generator
%   determined by the 'state' method, and initializes the state of that
%   generator using the value of STATE. RANDERR may not accept STATE in
%   a future release.  Use S instead.
%
%   This function can be useful for testing error-control coding.
%
%   Examples:
%
%   % Create a 2-by-5 error pattern with one error per row.
%   out = randerr(2,5)
%
%   % Create a 2-by-5 error pattern with two errors per row.
%   out = randerr(2,5,2)
%
%   % Create a 2-by-5 error pattern with one or three errors in each row.
%   out = randerr(2,5,[1 3])
%
%   % Create a 2-by-5 error pattern containing one error with 80% probability or
%   % three errors with 20% probability in each row.
%   out = randerr(2,5,[1 3;0.8 0.2])
%
%   See also RAND, RANDSRC, RANDI, RandStream.

%   Copyright 1996-2008 The MathWorks, Inc.
%   $Revision: 1.7.4.5 $  $Date: 2009/01/05 17:45:08 $


% Basic error checking and parameter setup etc.
error(nargchk(1,4,nargin,'struct'));

% --- Placeholder for the signature string.
sigStr = '';
n = [];
err = [];
state = [];

% --- Identify string and numeric arguments
isStream = false;
for i=1:nargin
   if(i>1)
      sigStr(size(sigStr,2)+1) = '/';
   end;
   % --- Assign the string and numeric flags
   if(isnumeric(varargin{i}))
      sigStr(size(sigStr,2)+1) = 'n';
   elseif(isa(varargin{i},'RandStream'))
      sigStr(size(sigStr,2)+1) = 'h';
      isStream = true;
   else
       error('comm:randerr:InvalidArg','Only numeric or RandStream arguments are accepted.');
   end;
end;

% --- Identify parameter signitures and assign values to variables
switch sigStr
   % --- randerr(m)
   case 'n'
      m		= varargin{1};

	% --- randerr(m, n)
	case 'n/n'
      m		= varargin{1};
      n		= varargin{2};

	% --- randerr(m, n, err)
	case 'n/n/n'
      m		= varargin{1};
      n  	= varargin{2};
      err	= varargin{3};

	% --- randerr(m, n, err, state)
	case {'n/n/n/n', 'n/n/n/h'}
      m		= varargin{1};
      n		= varargin{2};
      err	= varargin{3};
      state = varargin{4};

   % --- If the parameter list does not match one of these signatures.
   otherwise
      error('comm:randerr:InvalidSyntax','Syntax error.');
end;

if isempty(m)
   error('comm:randerr:EmptyM','Required parameter empty.');
end
if isempty(n)
   n = m;
end
if isempty(err)
   err = 1;
end

if (length(m) > 1) || (length(n) > 1)
   error('comm:randerr:NonScalarMatrixDims',...
       'Matrix dimensions must be scalars.');
elseif (floor(m) ~= m) || (floor(n) ~= n) || (~isreal(m)) || (~isreal(n))
   error('comm:randerr:NonIntMatrixDims',...
       'Matrix dimensions must be real integers.');
elseif (m <= 0) || (n <= 0)
   error('comm:randerr:NonPositiveMatrixDims',...
       'Matrix dimensions must be positive.');
elseif (~isfinite(m)) || (~isfinite(n))
   error('comm:randerr:InvalidInput',...
       'Input parameters must be finite.');
end

if any(~isfinite(err(1,:))) || any(~isreal(err(1,:))) ...
        || max(floor(err(1,:)) ~= err(1,:)) || max(err(1,:) < 0)
   error('comm:randerr:InvalidErrorsElements',...
       ['First row of the ERRORS parameter must contain finite real '...
       'positive integers.']);
elseif max(err(1,:) > n)
   error('comm:randerr:InvalidErrorsForm',...
       'Cannot place more "1"s in a row than there are elements in the row.');
end

[em, en] = size(err);

% If the probabilities are explicitly defined in the function call.
if (em == 2)

   if (~isreal(err(2,:)))
      error('comm:randerr:NonReal2ndRowErrorsElements',...
          'Second row of the ERRORS parameter must be real.');
   elseif any( (err(2,:) > 1) | (err(2,:) < 0) )
      error('comm:randerr:Invalid2ndRowErrorsVal',...
          ['Elements of second row of the ERRORS parameter must be '...
          'between zero and one.']);
   elseif ( abs(sum(err(2,:)) - 1) > sqrt(eps) )
      error('comm:randerr:InvalidProbabilitySum',...
          'Sum of ERRORS probability elements must equal one.');
   end

	% The 'prob' vector facilitates determining how many ones are to be
	% placed in each row of the output.
	% Example:
	%     If the 2nd row of 'err' is set to [0.1 0.2 0.3 0.4],
	%     then 'prob' will be [0.1 0.3 0.6 1.0].
   prob = cumsum(err(2,:));

% Default to using equal probabilities if not explicitly given.
elseif (em == 1)

   prob = (1:en) / en;

else
   error('comm:randerr:InvalidErrorsDims',...
       'ERRORS parameter cannot contain more than two rows.');
end

% Initialize the output matrix to all zeros.
out = zeros(m,n);

% Set the stream if specified.
if ~isempty(state)
    if ~isStream
        validateattributes(state, {'double', 'RandStream'}, ...
            {'real', 'scalar', 'integer'}, 'RANDERR', 'S');
        hStream = RandStream('mt19937ar', 'Seed', state);
    else
        hStream = state;
    end
    randFunc = @(a,b)rand(hStream, a, b);
else
    randFunc = @(a,b)rand(a, b);
end

% Now loop through each row of the output matrix.
for i = 1:m,

   % First determine how many ones will be put in each row.
   num = randFunc(1,1);
   not_done = 1;
   j = 1;
   while not_done
      if ( num < prob(j) )
         num = err(1,j);
         not_done = 0;
      end
      j = j + 1;
   end

   % Now find some random locations to place each one.
	[ignore,p] = sort(randFunc(1,n));
	out(i,:) = p <= num;

end;

% [EOF] randerr.m
