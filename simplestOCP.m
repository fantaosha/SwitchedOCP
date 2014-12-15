%SIMPLESTOCP Implements the most basic switched optimal control problem
%   The problem is a simple 1D problem with two modes, "up" or "down". The
%   running cost minimizes the distance to the origin.


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/


%% Data

% Relaxation order
order = 5; % Choose order >= 1;


%% Create problem definition structure
ocpDef.nModes = 2;
ocpDef.nStates = 1;
ocpDef.nControls = 0;
ocpDef.nLifts = 0;
ocpDef.scaling.t =1;
ocpDef.scaling.x =1;
ocpDef.scaling.u =1;
ocpDef.scaling.l =1;
ocpDef.dynamics{1} = @(t,x,u,l) 1;
ocpDef.dynamics{2} = @(t,x,u,l) -1;
ocpDef.runningCost{1} = @(t,x,u,l) x^2;
ocpDef.runningCost{2} = @(t,x,u,l) x^2;
ocpDef.initialCost = @(t,x,l) 0;
ocpDef.terminalCost = @(t,x,l) 0;
ocpDef.runningConstraints = @(t,x,u,l) [t*(1-t)>=0; 1-x^2>=0]; % time in [0,1], state in [-1 1]
ocpDef.initialConstraints = @(t,x,l) [x==0.5;t==0]; % start at x(t=0)=0.5
ocpDef.terminalConstraints = @(t,x,l) t==1; % do not forget to impose final time!
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


%% Solution extraction

% define extraction grid
npoints.t = 21;
npoints.var = 51; % je t'aime

% perform extraction on given grid with moment data as in the current
% Gloptipoly's state (viz., the last optimization performed)
[t,x,u,l,d] = extractSolution( ocpDef, measureSystem, order, npoints );

% plot reconstructed trajectory
figure
plot(t,x,'b*-');
xlabel('t');
ylabel('x');

% plot reconstructed modes
figure
plot(t,d(:,1),'bo-',t,d(:,2),'r+-');
xlabel('t');
ylabel('duty cycle');
legend('Mode 1','Mode 2');


