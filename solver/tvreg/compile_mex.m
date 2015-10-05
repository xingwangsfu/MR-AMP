function compile_mex
% This function compiles the MEX files included with tvreg.  Please
% configure the following compiler options for linking with FFTW.
%
% For help, there are several resources:
%    * The included PDF, section 7.3 on compiling instructions for MEX
%    * The FFTW website 
%          http://www.fftw.org/install/windows.html
%          http://www.fftw.org/install/install-Mac.html
%          http://www.fftw.org/install/install-Linux.html
%    * The MEX documentation
%          http://www.mathworks.com/support/tech-notes/1600/1605.html
%      Also try "help mex" on the MATLAB console.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configure linking with the FFTW3 library (http://www.fftw.org)
% If FFTW3 is not in MEX's search path, use the following variables.  Set
%    "libfftw3include"     to the location of fftw3.h, 
%    "libfftw3"            compiler options for linking with libfftw3
%    "libfftw3f"           compiler options for linking with libfftw3f
% libfftw3include = '"D:\libs\fftw"';
% libfftw3 = '"D:\libs\fftw\libfftw3-3.lib"';
% libfftw3f = '"D:\libs\fftw\libfftw3f-3.lib"';

% numsingle sets the floating point precision of the computations
%    numsingle = false       use double precision
%    numsingle = true        use single precision
numsingle = false;

% matlabctrlc enables an undocumented MEX method for Ctrl-C support.  If
% this is enabled, please also configure linking with the libut library.
%    matlabctrlc = false     disables Ctrl-C support
%    matlabctrlc = true      enables Ctrl-C support
matlabctrlc = true;

% Configure linking with the libut library
% Set the "libut" variable to the compiler options needed to link with the 
% libut library.  
%
% For a Windows machine, assuming the file is on MEX's search path
% libut = 'libut.lib';
% Otherwise, specify the full path, for example
% libut = ['"' matlabroot '\extern\lib\win32\microsoft\libut.lib"'];
% For a UNIX machine, this might work
% libut = '-lut';

% Set any compiler options to pass to MEX
opt = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code guesses the compiler options for anything that was not
% specified above.

if ~exist('libfftw3', 'var') || isempty(libfftw3) %#ok<NODEF>
    libfftw3 = '-lfftw3';
end

if ~exist('libfftw3f', 'var') || isempty(libfftw3f) %#ok<NODEF>
    libfftw3f = '-lfftw3f';
end

if ~exist('libut', 'var') || isempty(libut) %#ok<NODEF>
    if ispc  
        libut = 'libut.lib';
    else
        libut = '-lut';
    end
end

if ~exist('opt', 'var')
    opt = '';
end

% We now build the compiler options to pass to MEX

if matlabctrlc
    opt = [libut ' ' opt];
end

if ~numsingle
    opt = [libfftw3 ' ' opt];
else
    opt = [libfftw3f ' ' opt];
end

if exist('libfftw3include', 'var')
    opt = [sprintf('-I%s ', libfftw3include) opt];
end

if matlabctrlc
    opt = ['-DMATLAB_CTRL_C ' opt];
end

if numsingle
    opt = ['-DNUM_SINGLE ' opt];
end

files = {'tvreg.c', 'chanvese.c'};

fprintf(...
['\nThis function attempts to compile the MEX functions in the tvreg\n', ...
'package.  Please note, this script assumes that you have a C compiler\n', ...
'and that MATLAB has been configured to use it.\n\n', ...
'For help, please see the comments at the top of %s.m.\n'], ...
    mfilename);

for k = 1:length(files)
    mexstring = sprintf('mex %s %s', files{k}, opt);
    fprintf('\n\n%s\n\n', mexstring);

    success = true;
    % Attempt to compile (if it fails, "success" is set to false)
    eval(mexstring, 'success = false;');

    if ~success
        fprintf([ ...
'\nSorry, it seems there was a problem compiling %s.  Please edit\n', ...
'%s.m and correct the configuration options at the top of the file.\n\n'], ...
            files{k}, mfilename);
        break;
    end
end

if success
    fprintf('\nMEX functions compiled successfully.\n\n');
end
