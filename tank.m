%TANK Implements a switched tank filling example.
%   This is the tank problem as considered in
%      R. Vasudevan, H. Gonzalez, R. Bajcsy, S. S. Sastry. "Consistent
%      Approximations for the Optimal Control of Constrained Switched
%      Systems Part 2: An Implementable Algorithm". SIAM Journal on
%      Optimization and Control, 2013.
%   It is however attacked here by the moment technique developped in
%      M. Claeys, J. Daafouz, D. Henrion. "Modal occupation measures and
%      LMI relaxations for nonlinear switched systems control.". To be
%      published.


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/


%% Data

% Scaling factors: it is essential to scale variables so that everything
% falls within a unit box, otherwise moment problem will be poorly scaled
% (Consider moments of the Dirac measure located at 10...)
xscale = 5;
tscale = 10;
lscale = sqrt(5);

% Relaxation order
order = 3; % Choose order >= 1;

%% Create problem definition structure
ocpDef.nModes = 2;
ocpDef.nStates = 2; 
ocpDef.nControls = 0; 
ocpDef.nLifts = 2;
ocpDef.scaling.x = [xscale,xscale];
ocpDef.scaling.t = tscale;
ocpDef.scaling.u = 1;
ocpDef.scaling.l = [lscale,lscale];
ocpDef.dynamics{1} = @(t,x,u,l) [(1-lscale*l(1))*tscale/xscale;(lscale*l(1)-lscale*l(2))*tscale/xscale];
ocpDef.dynamics{2} = @(t,x,u,l) [(2-lscale*l(1))*tscale/xscale;(lscale*l(1)-lscale*l(2))*tscale/xscale];
ocpDef.runningCost{1} = @(t,x,u,l) 2*(x(2)*xscale-3)^2*tscale;
ocpDef.runningCost{2} = @(t,x,u,l) 2*(x(2)*xscale-3)^2*tscale;
ocpDef.initialCost = @(t,x,l) 0;
ocpDef.terminalCost = @(t,x,l) 0;
ocpDef.runningConstraints = @(t,x,u,l) [t*(1-t)>=0; % normalized time in [0,1]
                                        x(1) >= 0; % levels are positive
                                        x(2) >= 0;
                                        0==(lscale*l(1))^2-xscale*x(1); %algebraic lifts for square roots
                                        l(1)>=0;
                                        0==(lscale*l(2))^2-xscale*x(2);
                                        l(2)>=0;]; 
ocpDef.initialConstraints = @(t,x,l) [x==2/xscale;t==0];
ocpDef.terminalConstraints = @(t,x,l) [t==1];
% NB: to satisfy Putinar's theorem, there should be for each measure a ball
% constraint on all variables. Ignore this assumption at your own risks.
ocpDef.integralConstraints = {};

%% Construct and solve GloptiPoly problem
% (requires the installation of Yalmip and a SDP solver. We recommend Mosek
% for speed and SeDuMi for accuracy)

% Create GloptiPoly measure objects with default names (NB: resets all
% GloptiPoly states to zero). To impose specific names, construct measures
% by hand.
measureSystem = switchedMeasureSystem(ocpDef);

% Create GloptiPoly msdp object, modeling a given order moment relaxation
P = switchedRelaxation( ocpDef, measureSystem, order ); 

% Solve GloptiPoly problem
mset('yalmip',true); % default SDP solver of Yalmip will be called. If
%commented, Yalmip is not called and GloptiPoly looks for SeDuMi.
[status,obj,m] = msol(P);
obj


%% Solution extraction
npoints.t = 51;
npoints.var = 51;
% NB: extract moments of one order less, because higher order moments are
% very unprecise, since they are not constrained enough
[t,x,u,l,d] = extractSolution( ocpDef, measureSystem, order-1, npoints );

figure
plot(t,x,'*-');
xlabel('t');
ylabel('x');
legend('x_1','x_2');

figure
plot(t,d(:,1),'o-',t,d(:,2),'+-');
xlabel('t');
ylabel('duty cycle');
legend('Mode 1','Mode 2');

% create/update Bocop directory
toBocop('initTank',ocpDef,t,x,u,l,d);