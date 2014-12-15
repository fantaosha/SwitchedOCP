function checkOcpDef(ocpDef)
% Checks validity of optimal control problem definition
%
% The structure ocpDef must have the following fields with the required
% values:
%       -nModes: an integer with the number of modes (put 1 if this is a
%                regular OCP);
%       -nStates: an integer with the number of states; 
%       -nControls: an integer with the number of controls (put 0 if none
%                    are present);
%       -nLifts: an integer with the number of algebraic lifts (put 0 if
%                    none are present);
%       -scaling: a structure with the following fields:
%                   -x: a column vector with a scaling factor for each
%                       state
%                   -t: a scalar with a scaling factor for time
%                   -u: a column vector with a scaling factor for each
%                       control; put the scalar 1 if there is no control.
%                   -l: a column vector with a scaling factor for each
%                       algebraic lift; put the scalar 1 if there is no lift.
%       -dynamics: a cell array (one cell per mode) of function handles of
%                  respectively time, state, control and lift, that returns
%                  for each mode the scaled (in the unit box) vector field
%                  (as a column vector) associated to each mode.
%       -runningCost: a cell array (one cell per mode) of function handles
%                     of respectively time, state, control and lift, that
%                     returns for each mode the scaled (in the unit box)
%                     running cost (as a scalar) associated to each mode.
%                     Pu t the handle "@(t,x,u,l) 0" is there is no running
%                     cost for a given mode.
%        -initialCost: a function handle of respectively time, state, and
%                      lift, that returns the inital cost for starting in a
%                      given configuration. Put the handle "@(t,x,l) 0"
%                      if there is no initial cost.
%        -terminalCost: a function handle of respectively time, state, and
%                      lift, that returns the terminal cost for finishing
%                      in a given configuration. Put the handle "@(t,x,l) 0"
%                      if there is no terminal cost.
%        -runningConstraints: a function handle of respectively time, state,
%                             control and lift, that returns an array of
%                            (in)equalities that must be satisfied for each
%                            mode. The (in)equalities must be construted as
%                            would be done in GloptiPoly, assuming the
%                            inputs of the function handle are GloptiPoly
%                            mpol objects. We recommend adding ball
%                            constraints for each of the variables to
%                            satisfy Putinar's theorem.
%        -initialConstraints: a function handle of respectively time, state,
%                             and lift, that returns an array of
%                            (in)equalities that must be satisfied by the
%                             starting configuration.
%        -terminalConstraints: a function handle of respectively time, state,
%                             and lift, that returns an array of
%                            (in)equalities that must be satisfied by the
%                             starting configuration.
%        -integralConstraints: must be an empty object. Will hold in the
%                              future definition of integral constraints.
%
% If not, an error is raised. Note that the integrity of function
% definition is not checked. What is just checked in addition to correct
% structure definition are simple error checks on dimension.
%
% NB: all functions must be polynomials obviously!


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/

% first check structure of ocpDef
if ~isstruct(ocpDef)
    error('Input ocpDef must be a structure');
end

if ~all(isfield(ocpDef,{'nModes','nStates','nControls','nLifts','scaling',...
        'dynamics','runningCost','initialCost','terminalCost','runningConstraints',...
        'initialConstraints','terminalConstraints','integralConstraints'}))
    error('Structure ocpDef has missing fields');
end

if ~isstruct(ocpDef.scaling)
    error('ocpDef.scaling must be a structure');
end

if ~all(isfield(ocpDef.scaling,{'t','x','u','l'}))
    error('ocpDef.scaling has missing fields');
end

% Now simple error checks on dimension
if ~isnumeric(ocpDef.scaling.t) || ~all(isfinite(ocpDef.scaling.t(:))) || ...
        any(ocpDef.scaling.t(:)<=0) || numel(ocpDef.scaling.t)~=1
    ('ocpDef.scaling.t must be a strictly positive scalar');
end

if ~isnumeric(ocpDef.scaling.x) || ~all(isfinite(ocpDef.scaling.x(:))) || ...
        any(ocpDef.scaling.x(:)<=0) || size(ocpDef.scaling.x,2)~=1 ||size(ocpDef.scaling.x,2)~=ocpDef.nStates
    ('ocpDef.scaling.x must be a column vector of strictly positive values, one for each state');
end

if ocpDef.nControls > 0
    if ~isnumeric(ocpDef.scaling.u) || ~all(isfinite(ocpDef.scaling.u(:))) || ...
            any(ocpDef.scaling.u(:)<=0) || size(ocpDef.scaling.u,2)~=1 ||size(ocpDef.scaling.u,2)~=ocpDef.nControls
        ('ocpDef.scaling.u must be a column vector of strictly positive values, one for each control');
    end
else
    if ocpDef.scaling.u ~=1
        error('When no controls are defined, ocpDef.scaling.u must be 1');
    end
end

if ocpDef.nLifts> 0
    if ~isnumeric(ocpDef.scaling.l) || ~all(isfinite(ocpDef.scaling.l(:))) || ...
            any(ocpDef.scaling.l(:)<=0) || size(ocpDef.scaling.l,2)~=1 ||size(ocpDef.scaling.l,2)~=ocpDef.nLifts
        ('ocpDef.scaling.u must be a column vector of strictly positive values, one for each control');
    end
else
    if ocpDef.scaling.l ~=1
        error('When no controls are defined, ocpDef.scaling.u must be 1');
    end
end

if numel(ocpDef.dynamics)~=ocpDef.nModes
    error('Declared number of modes does not match size of ocpDef.dynamics');
end

if numel(ocpDef.runningCost)~=ocpDef.nModes
    error('Declared number of modes does not match size of ocpDef.runningCost');
end





