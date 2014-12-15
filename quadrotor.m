%QUADROTOR Implements a quadrotor example.
%   This is the quadrotor ("starmac") problem as considered in
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

M = 1.3;
g = 9.8;
L = 0.3050;
I = 0.0605;

% Scaling factors: it is essential to scale variables so that everything
% falls within a unit box, otherwise moment problem will be poorly scaled
% (Consider moments of the Dirac measure located at 10...)
xscale = [10; 2; 1; 1; 5; 0.2; 0.05]; % state order: [pos1,pos2,costheta,sintheta,pos1dot,pos2dot,thetadot]
tscale = 7.5;
uscale = 1e-3;

% Relaxation order
order = 2; % Choose order >= 2;

%% Create problem definition structure
ocpDef.nModes = 3;
ocpDef.nStates = 7; 
ocpDef.nControls = 1;
ocpDef.nLifts = 0;
ocpDef.scaling.x = xscale;
ocpDef.scaling.t = tscale;
ocpDef.scaling.u = uscale;
ocpDef.scaling.l = 1;
% state order: [pos1,pos2,costheta,sintheta,pos1dot,pos2dot,thetadot]
ocpDef.dynamics{1} = @(t,x,u,l) [xscale(5)*x(5)*tscale/xscale(1);
                                 xscale(6)*x(6)*tscale/xscale(2);
                                 -xscale(7)*x(7)*xscale(4)*x(4)*tscale/xscale(3);
                                 xscale(7)*x(7)*xscale(3)*x(3)*tscale/xscale(4);
                                 xscale(4)*x(4)/M*(uscale*u+M*g)*tscale/xscale(5);
                                 (xscale(3)*x(3)/M*(uscale*u+M*g)-g)*tscale/xscale(6); 
                                 0];
ocpDef.dynamics{2} = @(t,x,u,l) [xscale(5)*x(5)*tscale/xscale(1);
                                 xscale(6)*x(6)*tscale/xscale(2);
                                 -xscale(7)*x(7)*xscale(4)*x(4)*tscale/xscale(3);
                                 xscale(7)*x(7)*xscale(3)*x(3)*tscale/xscale(4);
                                 xscale(4)*x(4)*g*tscale/xscale(5);
                                 (xscale(3)*x(3)*g-g)*tscale/xscale(6);
                                 -L*uscale*u/I*tscale/xscale(7)];
ocpDef.dynamics{3} = @(t,x,u,l) [xscale(5)*x(5)*tscale/xscale(1);
                                 xscale(6)*x(6)*tscale/xscale(2);
                                 -xscale(7)*x(7)*xscale(4)*x(4)*tscale/xscale(3);
                                 xscale(7)*x(7)*xscale(3)*x(3)*tscale/xscale(4);
                                 xscale(4)*x(4)*g*tscale/xscale(5);
                                 (xscale(3)*x(3)*g-g)*tscale/xscale(6);
                                 L*uscale*u/I*tscale/xscale(7)];
running_cost = @(t,x,u,l) 5*(u*uscale)^2*tscale;
ocpDef.runningCost{1} = running_cost;
ocpDef.runningCost{2} = running_cost;
ocpDef.runningCost{3} = running_cost;
ocpDef.initialCost = @(t,x,l) 0;
ocpDef.terminalCost = @(t,x,l) (5*(xscale(1)*x(1)-6)^2 + 5*(xscale(2)*x(2)-1)^2 + (1-xscale(3)*x(3))/2); % NB (sin(theta/2))^2 = (1 - cos theta)/2
ocpDef.runningConstraints = @(t,x,u,l) [t*(1-t)>=0; % normalized time in [0,1]
                                        u*(1-u)>=0; % normalized control in [0,1]
                                        (xscale(3)*x(3))^2+(xscale(4)*x(4))^2==1; % attitude on SO2
                                        1-x.^2>=0]; % bound all reduced states on [-1,1]
ocpDef.initialConstraints = @(t,x,l) [t==0;
                                     x==[0;1;1;0;0;0;0]./xscale];
ocpDef.terminalConstraints = @(t,x,l) [t==1; % impose final scaled time
                                       1-x.^2>=0]; % bound all scaled states on [-1,1]
% NB: to satisfy Putinar's theorem, there should be for each measure a ball
% constraint on all variables and/or for each of them. Ignore this assumption at your own risks.
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
% NB: extract moments of one order less, because higher order moments are
% very unprecise, since they are not constrained enough
[t,x,u,l,d] = extractSolution( ocpDef, measureSystem, order, npoints );

% plot reconstructed trajectory
figure
plot(t,x,'*-');
xlabel('t');
ylabel('x');

% plot reconstructed control
figure
plot(t,u,'*-');
xlabel('t');
ylabel('u');

% plot reconstructed modal duty cycles
figure
plot(t,d(:,1),'o-',t,d(:,2),'+-',t,d(:,3),'*-');
xlabel('t');
ylabel('duty cycle');
legend('Mode 1','Mode 2','Mode3');

% create/update Bocop directory 
toBocop('initQuadrotor',ocpDef,t,x,u,l,d);
