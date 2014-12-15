function toBocop(bocopdir,ocpDef,t,x,u,l,d)
%TOBOCOP Create/update folder with BOCOP initialization files 
% toBocop(bocopdir,ocpDef,t,x,u,d) creates/updates folder bocopdir with
%   BOCOP initialisation files for optimal control problem described by
%   ocpDef. Folder bocopdir must be the "init" folder of the corresponding
%   BOCOP implementation of the problem described in ocpDef.
%
%   Time series for time stamps, states, control, lifts and duty
%   cycles are given by resp. t, x, u, l and d. Each line of these arrays
%   must correspond to the same time stamp. The different init files
%   created follow the conventation:
%       - the n states defined in ocpDef are printed in the state.[i].init
%       file, where [i] ranges from 0 to n-1 (beware of that c++ indexing).
%       - the state.[n+1].init file is an approximation of the running
%       cost, which must be entered as an additional state in BOCOP.
%       - algebraic lifts are not printed
%       - the controls, if any, are printed in the control.[i].init files
%       - modes are modelled as additional controls and appended to the
%       control.[i].init files. In BOCOP, they must be constrained to add
%       up to 1 for every time.
%
%   Note that currently, there is a bug in BOCOP's handling of input times,
%   which must be scaled in the [0,1] interval. The function automatically
%   scales time to that interval by dividing t by its largest value


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/


%% Input parsing

% create folderif needed
 if exist(bocopdir,'file')~=7
    mkdir(bocopdir);
 end
 
% check consistency of inputs

% check ocpDef
checkOcpDef(ocpDef);

% check t
if ~isnumeric(t) || ~all(isfinite(t(:))) || size(t,2)~=1
    error('t must be a column vector of finite numeric vales');
end

% check x
if ~isnumeric(x) || ~all(isfinite(x(:))) || size(x,1)~=length(t) || size(x,2)~=ocpDef.nStates
    error('x must be a matrix of finite numeric vales, with as many lines as t and a column for each state.');
end

% check u
if ocpDef.nControls > 0
    if ~isnumeric(u) || ~all(isfinite(u(:))) || size(u,1)~=length(t) || size(u,2)~=ocpDef.nControls
        error('u must be a matrix of finite numeric vales, with as many lines as t and a column for each contol.');
    end
end

% check l
if ocpDef.nLifts > 0
    if ~isnumeric(l) || ~all(isfinite(l(:))) || size(l,1)~=length(t) || size(l,2)~=ocpDef.nLifts
        error('l must be a matrix of finite numeric vales, with as many lines as t and a column for each lift.');
    end
end

% check d
if ~isnumeric(d) || ~all(isfinite(d(:))) || size(d,1)~=length(t) || size(d,2)~=ocpDef.nModes
    error('d must be a matrix of finite numeric vales, with as many lines as t and a column for each mode.');
end


% compute running cost
% store here for each time the integrand
dcost = zeros(size(t)); 
% compute the integrand value
for k=1:length(t)
    for im=1:ocpDef.nModes
        dcost(k) = dcost(k) + ...
            d(k,im)*ocpDef.runningCost{im}(t(k)/ocpDef.scaling.t,x(k,:)./ocpDef.scaling.x(:)',u(k,:)./ocpDef.scaling.u(:)',l(k,:)./ocpDef.scaling.l(:)');
    end
end
% use simple trapezoidal rule
runningCost = cumtrapz( t/ocpDef.scaling.t, dcost );

%% print all necessary files

%%% TODO remove scaling hack to counter BOCOP's bug
tscale = max(t);

% states
for is = 1:ocpDef.nStates
    filename = fullfile(bocopdir,sprintf('state.%d.init',is-1)); % -1 because Bocop indexing is zero-based
    initBocopFile( filename, t/tscale, x(:,is));
    dlmwrite(fullfile(bocopdir,sprintf('state.%d.txt',is-1)), [t x(:,is)]);
end

% running cost
filename = fullfile(bocopdir,sprintf('state.%d.init',ocpDef.nStates));
initBocopFile( filename, t/tscale, runningCost);
dlmwrite(fullfile(bocopdir,sprintf('state.%d.txt',ocpDef.nStates)), [t runningCost]);

% controls
for ic = 1:ocpDef.nControls
    filename = fullfile(bocopdir,sprintf('control.%d.init',ic-1)); % -1 because Bocop indexing is zero-based
    initBocopFile( filename, t/tscale, u(:,ic));
    dlmwrite(fullfile(bocopdir,sprintf('control.%d.txt',ic-1)), [t u(:,ic)]);
end

% duty cycles
for im = 1:ocpDef.nModes
    filename = fullfile(bocopdir,sprintf('control.%d.init',ocpDef.nControls-1+im)); % -1 because Bocop indexing is zero-based
    initBocopFile( filename, t/tscale, d(:,im));
    dlmwrite(fullfile(bocopdir,sprintf('control.%d.txt',ocpDef.nControls-1+im)), [t d(:,im)]);
end


end

function initBocopFile( filename, t, x, varname)
%INITBOCOPFILE Create BOCOP init file from time series
% initBocopFile( filename, t, x), where filename is a string
%   containing a target filename (with possible directory path), t is an
%   array of time stamps and x the corresponding data value, creates an
%   init file compatible with BOCOP (www.BOCOP.org). Time series do not
%   need to be ordered or be in vector form.
% initBocopFile( filename, t, x, varname) uses given variable name in the
%   file's header. Default value is deduced from filename using Matlab's
%   fileparts.
%
% Example
%   initBocop( 'state.0.init', 0:.1:10, sin(0:.1:10))


%% error checking and input parsing
if ~isfloat(t) || ~isreal(t) || ~all(isreal(t)) || isempty(t)
   error('Input t must be a valid numeric array');
end

if ~isfloat(x) || ~isreal(x) || ~all(isreal(x)) || isempty(x)
   error('Input x must be a valid numeric array');
end

if numel(t)~=numel(x)
    error('Inputs t and x must have the same number of elements');
end

% reshape and sort arrays
t = t(:);
x = x(:);
[t,idx] = sort(t);
x = x(idx);
np = length(t);

if ~ischar(filename)
    error('File name must be a string'); 
end

if nargin <= 3
    [dummy,varname] = fileparts(filename);
else
    if ~ischar(varname)
       error('Variable name must be a string'); 
    end
end

fid = fopen(filename,'w');
if fid == -1
   error('Error opening file'); 
end


%% Print array in file
fprintf(fid,'# Starting point file.\n');
fprintf(fid,'# This file contains the values of the initial points\n');
fprintf(fid,'# for variable %s\n',varname);
fprintf(fid,'\n');

fprintf(fid,'# Type of initialization :\n');
fprintf(fid,'linear\n');
fprintf(fid,'\n');
fprintf(fid,'# Number of interpolation points :\n');
fprintf(fid,'%d\n',np);
fprintf(fid,'\n');
fprintf(fid,'# Interpolation points for the starting point :\n');
for k = 1:np
    fprintf(fid,'%.16f %.16f\n',t(k),x(k));
end

fclose(fid);

end

