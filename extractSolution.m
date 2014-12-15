function [tgrid,x,u,l,d,weights] = extractSolution( ocpDef, measDef, order, npoints)
%EXTRACTSOLUTION Creates folder with BOCOP initialization files 
% [t,x,u,d] = extractSolution( ocpDef, measDef, order, npoints )
%   Extract solution (column vector t of time stamps, nstamps x nstates
%   array x with corresponding state values, u for controls, d for modal
%   duty cycles) from the moments of measure system measDef up to given
%   order, with given number of points (given in as tructure with fields t
%   and var).
%
%   Note that modal duty cycles sum up to 1 for (almost) all time by
%   definition.
%
%   Note that it is implied that the problem is a fixed time OCP, scaled so
%   that normalized time falls into the [0,1] interval and normalized
%   controls and states fall each into [-1,1]. The returned data, however,
%   is not normalized.
%
%   The function also makes the crucial assumption that the optimal
%   solution is unique. If not, the problem must be approprately perturbed.
%
%   The function requires an LP solver called linprog with the same calling
%   syntax as the Optimization Toolbox' linprog.
%
%   See Mathieu Claeys and Rodolphe Sepulchre, "Reconstructing trajectories
%   from the moments of occupation measures", CDC '14.


% Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/

% generate grid
tgrid = linspace(0,1,npoints.t)';
xgrid = linspace(-1,1,npoints.var)';

x = zeros(length(tgrid),ocpDef.nStates);
u = zeros(length(tgrid),ocpDef.nControls);
l = zeros(length(tgrid),ocpDef.nLifts);
d = zeros(length(tgrid),ocpDef.nModes);

% loop over all states
for is=1:ocpDef.nStates
    
    % extract weight matrix
    W = identifyVariable( measDef.modal.x(is,:), measDef, tgrid, xgrid, order);
    weights.x(is).W = W;
    
    % average for all time
    for k = 1:length(tgrid)-1
        for im = 1:ocpDef.nModes
            x(k,is) = x(k,is)+xgrid'*W{im}(:,k)*ocpDef.scaling.x(is);
        end
    end
    % repeat final value, since fitting is done on intervals and not at
    % nodes
    x(end,is) = x(end-1,is);
     
end


% loop over all controls
for ic=1:ocpDef.nControls
    
    % extract weight matrix
    W = identifyVariable( measDef.modal.u(ic,:), measDef, tgrid, xgrid, order);
    weights.u(im).W = W;
    
    % average for all time
    for k = 1:length(tgrid)-1
        for im = 1:ocpDef.nModes
            u(k,ic) = u(k,ic)+xgrid'*W{im}(:,k)*ocpDef.scaling.u(ic);
        end
    end
    % repeat final value, since fitting is done on intervals and not at
    % nodes
    u(end,ic) = u(end-1,ic);
     
end


% loop over all lifts
for il=1:ocpDef.nLifts
    
    % extract weight matrix
    W = identifyVariable( measDef.modal.l(il,:), measDef, tgrid, xgrid, order);
    
    % average for all time
    for k = 1:length(tgrid)-1
        for im = 1:ocpDef.nModes
            l(k,il) = l(k,il)+xgrid'*W{im}(:,k)*ocpDef.scaling.l(il);
        end
    end
    % repeat final value, since fitting is done on intervals and not at
    % nodes
    l(end,il) = l(end-1,il);
     
end


% treat modal duty cycles
d = extractDutyCycles( measDef, tgrid, order );
% repeat final value
d = [d;d(end,:)];


% rescale time
tgrid = tgrid*ocpDef.scaling.t;

end


function [W,T,X] = identifyVariable( theVars, measDef, tgrid, vargrid, order)


[dummy,pow] = genind(2,2*order);
expT = pow(:,1);
expX = pow(:,2);


nm = length(measDef.modal.measure);
ny = length(expT);
nt = numel(tgrid)-1;
nx = numel(vargrid);
nc = nt*nx;

midpoints = tgrid(1:end-1)+diff(tgrid)/2;
[T,X] = meshgrid(midpoints,vargrid);

%% Solve LP


% define moments of cells
momcells = zeros(ny,nc);

for iy = 1:ny
    for ic = 1:nc
        [idxX,idxT] = ind2sub([nx nt],ic); % order of outputs matches meshgrid's numbering convention
        momcells(iy,ic) = (tgrid(idxT+1)^(expT(iy)+1)-tgrid(idxT)^(expT(iy)+1))/(expT(iy)+1) * vargrid(idxX)^expX(iy);
    end
end

A = zeros(2*nm*ny,nc*nm+1);
b = zeros(2*nm*ny,1);
for j=1:nm
    A(ny*(j-1)+(1:ny),nc*(j-1)+(1:nc)) = momcells;
    A(ny*(j-1)+(1:ny),end) = -1;
    A(ny*nm+ny*(j-1)+(1:ny),nc*(j-1)+(1:nc)) = -momcells;
    A(ny*nm+ny*(j-1)+(1:ny),end) = -1;
    b(ny*(j-1)+(1:ny)) = double(mom(measDef.modal.t(j).^expT.*theVars(j).^expX));
    b(ny*nm+ny*(j-1)+(1:ny)) = -b(ny*(j-1)+(1:ny));
end


Aeqbloc = zeros(nt,nc);
for it=1:nt
    Aeqbloc(it,(it-1)*nx+(1:nx)) = 1;
end

Aeq = [repmat(Aeqbloc,1,nm) zeros(nt,1)];
beq = ones(nt,1);

c= zeros(1,nc*nm+1);
c(end) = 1;

[sol, fval] = linprog(c,A,b,Aeq,beq,zeros(nc*nm+1,1));


W=cell(nm,1);
for im=1:nm
    W{im} = reshape ( sol((im-1)*nc+(1:nc)), nx, nt );
end



end



function [W,fval] = extractDutyCycles( measDef, tgrid, order )
%EXTRACTDUTYCYCLE Extract duty cycles of different modes from moment data
% W = extractDutyCycles( measDef, T, order ), where
%   measDef is a structure defining a switched system GloptiPoly measure
%   system, T is a grid of time stamps in ascending order and order is the
%   given extraction order, returns a n_modes x n_timestamps-1 matrix W of
%   duty cycle fitting in an optimal way the observed numeric moments of
%   the measure system. Moments up to the given order must have been
%   previously computed via GloptiPoly. W is scalled so that for each time
%   interval, the corresponding duty cycles add up to one.
%
% [W,normfit] = extractDutyCycles(...) also returns the norm of the fitting
%   error.
%
%   The function requires a LP solver called linprog with the same calling
%   syntax as the Optimization Toolbox' linprog.


%% Error checking and input parsing

if exist('linprog','file')~=2
    error('Linear programming routine required');
end

if ~isfloat(tgrid) || ~isreal(tgrid) || ~all(isreal(tgrid)) || isempty(tgrid)
   error('Input t must be a valid numeric array');
end
if numel(tgrid) < 2
   error('Time grid must have at least 2 elements');
end

nm = length(measDef.modal.measure);
ny = 2*order+1;
ne = numel(tgrid)-1;

exponents = (0:2*order)';

%% Construct LP
Y = zeros(ny,ne);

for ie = 1:length(exponents)
    for ife = 1:ne
        Y(ie,ife) = tgrid(ife+1)^(exponents(ie)+1)/(exponents(ie)+1) - tgrid(ife)^(exponents(ie)+1)/(exponents(ie)+1);
    end
end

A = zeros(2*ny*nm,ne*nm+1);
b = zeros(2*ny*nm,1);
for j=1:nm
    A(ny*(j-1)+(1:ny),ne*(j-1)+(1:ne)) = Y;
    A(ny*(j-1)+(1:ny),end) = -1;
    A(ny*nm+ny*(j-1)+(1:ny),ne*(j-1)+(1:ne)) = -Y;
    A(ny*nm+ny*(j-1)+(1:ny),end) = -1;
    b(ny*(j-1)+(1:ny)) = double(mom(measDef.modal.t(j).^exponents));
    b(ny*nm+ny*(j-1)+(1:ny)) = -double(mom(measDef.modal.t(j).^exponents));
end


Aeq = [repmat(eye(ne),1,nm) zeros(ne,1)];
beq = ones(ne,1);

c= zeros(1,ne*nm+1);
c(end) = 1;

[sol, fval] = linprog(c,A,b,Aeq,beq,zeros(ne*nm+1,1));



W = reshape( sol(1:end-1),  ne, nm );




end


