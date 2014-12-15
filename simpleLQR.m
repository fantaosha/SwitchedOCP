%SIMPLELQR Implements a simple switched LQR problem.
%   This is the LQR problem as considered in
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

% For dynamics definition
A = [1.0979 -0.0105 0.0167; -0.0105 1.0481 0.0825; 0.0167 0.0825 1.1540];
B1 = [0.9801; -0.1987; 0];
B2 = [0.1743; 0.8601; -0.4794];
B3 = [0.0952; 0.4699; 0.8776];

% Scaling factors: it is essential to scale variables so that everything
% falls within a unit box, otherwise moment problem will be poorly scaled
% (Consider moments of the Dirac measure located at 10...)
xscale = 2;
tscale = 2;
uscale = 20;

% Relaxation order
order = 2; % Choose order >= 1;

%% Create problem definition structure
ocpDef.nModes = 3;
ocpDef.nStates = 3;
ocpDef.nControls = 1;
ocpDef.nLifts = 0;
ocpDef.scaling.x = [xscale,xscale,xscale];
ocpDef.scaling.t = tscale;
ocpDef.scaling.u = uscale;
ocpDef.scaling.l = 1;
ocpDef.dynamics{1} = @(t,x,u,l) A*x*tscale+B1*u*uscale*tscale/xscale;
ocpDef.dynamics{2} = @(t,x,u,l) A*x*tscale+B2*u*uscale*tscale/xscale;
ocpDef.dynamics{3} = @(t,x,u,l) A*x*tscale+B3*u*uscale*tscale/xscale;
ocpDef.runningCost{1} = @(t,x,u,l) 0.01*(u*uscale)^2*tscale;
ocpDef.runningCost{2} = @(t,x,u,l) 0.01*(u*uscale)^2*tscale;
ocpDef.runningCost{3} = @(t,x,u,l) 0.01*(u*uscale)^2*tscale;
ocpDef.initialCost = @(t,x,l) 0;
ocpDef.terminalCost = @(t,x,l) sum((xscale*x-1).^2);
ocpDef.runningConstraints = @(t,x,u,l) [0<=t*(1-t); % normalized time in [0,1]
                                        0<=1-x.^2;
                                        0<=1-u.^2]; 
ocpDef.initialConstraints = @(t,x,l) [0==x;0==t];
ocpDef.terminalConstraints = @(t,x,l) [1==t; 0<=1-x.^2];
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

% define extraction grid
npoints.t = 51;
npoints.var = 51;

% perform extraction on given grid with moment data as in the current
% Gloptipoly's state (viz., the last optimization performed)
[t,x,u,l,d] = extractSolution( ocpDef, measureSystem, order, npoints);

% plot reconstructed trajectory
figure
plot(t,x,'*-');
xlabel('t');
ylabel('x');
legend('x_1','x_2','x_3');

% plot reconstructed duty cycle
figure
plot(t,d(:,1),'o-',t,d(:,2),'+-',t,d(:,3),'*-');
xlabel('t');
ylabel('duty cycle');
legend('Mode 1','Mode 2','Mode 3');

% create/update Bocop directory
toBocop('initLQR',ocpDef,t,x,u,l,d);
