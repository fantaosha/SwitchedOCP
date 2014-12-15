This folder contains the code to run the examples presented in

M. Claeys, J. Daafouz, D. Henrion. "Modal occupation measures and
LMI relaxations for nonlinear switched systems control". To be
published.


------------------- Copyright Notice ----------------------------
Copyright 2014 Mathieu Claeys, http://mathclaeys.wordpress.com/

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


------------------- Folder composition ----------------------------
A bare-bone switched optimal control toolbox composed of the following
functions:
   -switchedmeasureSystem.m: creates all necessary GloptiPoly objects
   -switchedRelaxation.m: creates the SDP relaxation
   -checkOcpDef.m: checks integrity of problem's definition
   -extractSolution.m: extracts approximate global candidate
   -toBocop.m: creates a folder of initalization file for Bocop

A set of working examples to base your own problems on, implemented as
scripts:
   -simplestOCP.m: implements the most basic switched OCP
   -simpleLQR.m: implements the LQR problem as in the given reference
   -tank.m: implements the double tank problem as in the given reference
   -quadrotor.m: implements the quadrotor problem as in the given reference

To run one of the examples, follow the installation instructions, then clear
the Matlab workspace and run the example script of your choosing. To
select relaxation order and SDP solver, you will need to
modify the scripts.

------------------- Installation ----------------------------
The following are required to run the scripts:
    - A recent Matlab install
    - GloptiPoly 3, a toolbox for constructing moment relaxation, which
      you can download at
                      http://homepages.laas.fr/henrion/software/gloptipoly/
    - A SDP solver (see below)

By default, GloptiPoly calls the SDP solver SeDuMi, which you can find at
                      http://sedumi.ie.lehigh.edu/ 
However, the user may want to solve the relaxations with a different
solver than SeDuMi, given the wide range of numerical stability and
speed between different solvers. To run other SDP solvers, GloptiPoly
uses Yalmip as an interface between moment relaxations and SDP solvers.
You will need then to download Yalmip at
                      http://users.isy.liu.se/johanl/yalmip/
You will be thus able to solve the moment relaxations via any SDP
solver that Yalmip can handle. We warmly recommand Mosek as the default
SDP solver, which you can find at
                      http://www.mosek.com/

To sum up, if you are happy with SeDuMi's speed then you only need the
basic install:
    -Matlab
    -GloptiPoly
    -SeDuMi
If you want to use Mosek (or any other SDP solver), then you need to
install:
    -Matlab
    -GloptiPoly
    -Yalmip
    -Mosek (or any other SDP solver)
Note that the example scripts suppose that you have installed the second
configuration. To revert to the first one, you need to comment out
the appropriate line which is clearly marked (or ctrl+F: "mosek").


-------------------- Using the toolbox -----------------------
The example scripts are pretty much self-explanatory, which should
allow you to build easily your own examples. Basically, you define
your problem in an appropriate structure, and the shipped function
allows you to automate all the tedious GloptiPoly code. Obviously,
implementing the GloptiPoly problem "by hand" will lead to faster
code.

The toolbox allows you to define "lifts". Those are necessary when your
problem is not polynomial per se, but can lifted as such, e.g.
    -sqrt(x), as treated in the double tank example
    -abs(x), which you substitute by l, and constrain by l^2=x^2 and l>=0
    -fraction x/y which you substitute by l, and constrain as ly=x
    -min, max, and many other examples (really!).
Obviously, lifts add variables to each measure, so they must be counted
as additional states in the complexity estimates of the method. There is
no such thing as a free lunch.

All scripts are followed by an example extraction procedure, to reconstruct
optimal trajectories, controls, and modal duty cycles. The procedure
requires a linprog routine (as implemented in Matlab's Optimization Toolbox)
or in Mosek's Matlab interface.

Those reconstructed objects are then plotted for inspection, and converted
into BOCOP init file. BOCOP is a free generic optimal control solver that
we warmly recommend, and which you can find at
                 http://bocop.saclay.inria.fr/
The scripts generate a folder "initNameOfProblem", with files ready to be
used by BOCOP as starting guess. Those files must be placed in the
"init" folder of the corresponding BOCOP's implementation of the problem.
We haven't implemented an automatic translator yet, so if you wish to try
this functionality, you will need to re-code the problem in BOCOP (fairly
easy though), following the conventions detailed if you type "help toBocop"
in Matlab.

WARNING: during development of these scripts, a bug was found in BOCOP for
the handling of time stamps of the init files. Newer versions of BOCOP
will correct this, and the function toBocop will need to be updated.

--------------Tips for implementing your own examples-------------------
We cannot stress enough two crucial points to "make it work":
    -If you are unsure of what you are doing, always constrain all
     variables of each measure quadratically, either individually or
     in a big sphere constraints, e.g:
               1-x(1)^2 >= 0 % constraints x(1) in [-1,1] interval
               u(2)*(1-u(2)) >= 0 % constraints u(2) in [0,1] interval
               1-sum(x.^2) >= 0 % all state in unique ball costraint
     The reason is to comply with the conditions of Putinar's theorem.
     Also, constraints add marginal computation times to the problem,
     but allow for richer (read: less conservative...) positivity
     certificates (when "algebraically independent", of course...).
    -Scale your problem so that all trajectories, controls, etc. fall
     within a unit box. Otherwise, the moment relaxations will be
     so poorly conditioned that only the lowest relaxations will
     suceed.

Also, check the following before calling for help:
    -The toolbox does not impose initial time constraint. This must be
     specified as an initialConstraint, e.g. t==0 
    -Same thing for the final time, e.g. t==1
    -Same thing for the initial conditions, e.g. x(1)==0 and x(2)==1
    -Problems with a state space above 6 become challenging, even for
     Mosek. There is nothing we can do about it.
    -It is very easy to run into "out of memory" Matlab errors when
     you have a 32 bit machine. This limitation is induced by Matlab
     and there is nothing we can do about it. You will need to borrow
     someone else's 64 bit machine to work around this problem.
    -Mosek is really much faster than SeDuMi...
    -Remember that GloptiPoly was designed to solve finite-dimensional
     optimization problems. The messages printed on the command window
     might then be meaningless in the context of optimal control. Also,
     the solution extraction that is performed by GloptiPoly is useless
     in our case, and you may want to modify the code of GloptiPoly to
     avoid performing it (sometimes, the code gets stucked there)