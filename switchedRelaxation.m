function P = switchedRelaxation( ocpDef, measDef, order )
%SWITCHEDRELAXATION Build moment relaxation given switched OCP
% P = switchedRelaxation( ocpDef, measDef, order ) build a GloptiPoly's
%    msdp object P containing the relaxation of order order of switched
%    optimal control problem ocpDef using measure system measDef. See
%    switchedMeasureSystem for a description of the various objects
%
% Note that the code is written so that no usefull moment substitution is
% induced. This is simply to avoid those "cyclic substitution" GloptiPoly
% errors. Again, to be more efficient, you will need to code that yourself,
% or modify carefully this function.


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/


%Define cost to minimize
argcrit = mom(   ocpDef.initialCost(measDef.initial.t, measDef.initial.x, measDef.initial.l),        indmeas(measDef.initial.measure)) + ...
          mom(   ocpDef.terminalCost(measDef.terminal.t, measDef.terminal.x, measDef.terminal.l),    indmeas(measDef.terminal.measure));
for j=1:ocpDef.nModes % loop over all modal measures
    argcrit = argcrit + mom(   ocpDef.runningCost{j}(measDef.modal.t(:,j), measDef.modal.x(:,j), measDef.modal.u(:,j), measDef.modal.l(:,j)),    indmeas(measDef.modal.measure(j))); 
end


%Define test functions of time and state for weak dynamics moment constraint
% first define full order test functions, which is of max degree 2*order
vi = mmon( [measDef.initial.t; measDef.initial.x], 2*order );
vt = mmon( [measDef.terminal.t; measDef.terminal.x], 2*order );
vm = cell([ocpDef.nModes,1]);
dvm = cell(size(vm));
for j=1:ocpDef.nModes % loop over all modal measures
    vm{j} = mmon( [measDef.modal.t(:,j); measDef.modal.x(:,j)], 2*order );
    dvm{j} = diff(vm{j},measDef.modal.t(:,j)) + ...
        diff(vm{j},measDef.modal.x(:,j))*ocpDef.dynamics{j}(measDef.modal.t(:,j), measDef.modal.x(:,j), measDef.modal.u(:,j), measDef.modal.l(:,j));
end
% then prune any constraint of too high of a degree
toRemove = false(size(vi)); % mask of test functions to remove because they exceed the required order
for k = 1:numel(toRemove) %loop over all test functions
    for j=1:ocpDef.nModes % loop over all modal measures
        if deg(dvm{j}(k)) > 2*order
            toRemove(k) = true;
            break; % no need to loop over further measures
        end
    end
end
vi(toRemove) = [];
vt(toRemove) = [];
for j=1:ocpDef.nModes % loop over all modal measures
    vm{j}(toRemove) = [];
    dvm{j}(toRemove) = [];
end


%Define moment constraints
momconstr = {};
% mass constraints
momconstr{end+1} = 1 == mass(measDef.initial.measure);
momconstr{end+1} = 1 == mass(measDef.terminal.measure);
% weak dynamics
summomsmodal = 0;
for j=1:ocpDef.nModes
    summomsmodal = summomsmodal + mom(dvm{j});
end
momconstr{end+1} = 0  == mom(vt)-mom(vi)-summomsmodal; %put everything to the right to avoid cycling problems
                                                       % For specific problems, this may be optimized
                                                       
% integral constraints: experimental only, use at own risk
if ~isempty(ocpDef.integralConstraints)
    % do some error checking
    if size(ocpDef.integralConstraints.Aeq,2)~=ocpDef.nModes
        error('Aeq must have as many columns as modes');
    end
    if size(ocpDef.integralConstraints.Aeq,1)~=numel(ocpDef.integralConstraints.beq)
        error('Aeq must have as many rows as elements in beq');
    end
    for ic=1:numel(ocpDef.integralConstraints.beq)
        rightmember = 0;
        for im=1:ocpDef.nModes
            rightmember =  rightmember + mom(ocpDef.integralConstraints.Aeq{ic,im}(measDef.modal.t(:,im), measDef.modal.x(:,im), measDef.modal.u(:,im), measDef.modal.l(:,im)));
        end
        momconstr{end+1} = ocpDef.integralConstraints.beq(ic)==rightmember;
    end
end
                                                       
                                                       
%Define support constraints
suppconstr = {};
arg = ocpDef.initialConstraints(measDef.initial.t, measDef.initial.x, measDef.initial.l);
if ~isempty(arg)
    suppconstr{end+1} = arg;
end
arg = ocpDef.terminalConstraints(measDef.terminal.t, measDef.terminal.x, measDef.terminal.l);
if ~isempty(arg)
    suppconstr{end+1} = arg;
end
for j=1:ocpDef.nModes
    arg = ocpDef.runningConstraints(measDef.modal.t(:,j), measDef.modal.x(:,j), measDef.modal.u(:,j), measDef.modal.l(:,j));
    if ~isempty(arg)
        suppconstr{end+1} = arg;
    end
end


% assemble problem
P = msdp(min(argcrit),momconstr{:},suppconstr{:},order);


end

