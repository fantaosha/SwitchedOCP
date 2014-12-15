function measDef = switchedMeasureSystem( ocpDef )
%SWITCHEDMEASURESYSTEM Get measure system from switched OCP description
% measDef = switchedMeasureSystem( ocpDef ) returns a structure pointing to
%    various GloptiPoly objects describing switched system optimal control
%    problems. See ocpdefCheck and the various shipped examples for the
%    format of ocpDef.
%
%    Note that the function resets the GloptiPoly internal state.
%
%See also ocpDefCheck, simpleLQR, tank, quadrotor


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/


%% Setup
checkOcpDef(ocpDef)

if exist('mset','file')~=2
    error('GloptiPoly is not installed properly');
end

%% Reset GloptiPoly state
mset('clear');

%% Create variables

% initial measure
mpol('ti',1);
measDef.initial.t = ti;
mpol('xi',ocpDef.nStates);
measDef.initial.x = xi;
if ocpDef.nLifts >= 1
    mpol('li',ocpDef.nStates);
    measDef.initial.l = li;
else
    measDef.initial.l = [];
end
measDef.initial.measure = meas(1);

% final measure
mpol('tt',1);
measDef.terminal.t = tt;
mpol('xt',ocpDef.nStates);
measDef.terminal.x = xt;
if ocpDef.nLifts >= 1
    mpol('lf',ocpDef.nStates);
    measDef.terminal.l = lf;   
else
    measDef.terminal.l = [];
end
meas([measDef.terminal.t;measDef.terminal.x;measDef.terminal.l]);
measDef.terminal.measure = meas(2);

% modal measure
mpol('tm',1,ocpDef.nModes);
measDef.modal.t = tm;
mpol('xm',ocpDef.nStates,ocpDef.nModes);
measDef.modal.x = xm;
if ocpDef.nControls >= 1
    mpol('um',ocpDef.nControls,ocpDef.nModes);
    measDef.modal.u = um;
else
    measDef.modal.u = zeros(0,ocpDef.nModes);
end
if ocpDef.nLifts >= 1
    mpol('lm',ocpDef.nLifts,ocpDef.nModes);
    measDef.modal.l = lm;
else
    measDef.modal.l = zeros(0,ocpDef.nModes);
end

for j=1:ocpDef.nModes
    meas([measDef.modal.t(:,j);measDef.modal.x(:,j);measDef.modal.u(:,j);measDef.modal.l(:,j)]);
    measDef.modal.measure(j) = meas(2+j);
end


end

