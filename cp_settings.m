%% Settings for the cart-pole scenario
%
% Based on a script by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
%% High-Level Steps
% 1 Define state and important indices
% 2 Set up scenario
% 3 Set up the plant structure
% 4 Set up the policy structure
% 5 Set up the cost structure
% 6 Set up the GP dynamics model structure
% 7 Parameters for policy optimization
% 7 Plotting verbosity
% 9 Some array initializations

%% 0. Preparation [optional]

% init number generators to make the experiments reproducible
rand('seed',1); 
randn('seed',1); 

% printing
format short; 
format compact; 

% include some paths
try
  rd = 'C:\Users\antal\sztaki\MATLAB\pilcoV0.9\';
  addpath([rd 'base'],[rd 'util'],[rd 'gp'],[rd 'control'],[rd 'loss']);
catch
end
basename = 'Cart_Pole_';          % filename used for saving data

%% 1. Define state and important indices

% 1a. State vector (including all augmentations)
%
%  1  x          cart position
%  2  v          cart velocity
%  3  dtheta     angular velocity of the pendulum
%  4  theta      angle of the pendulum
%  5  sin(theta) cartesian representation of theta [augmentation]
%  6  cos(theta) cartesian representation of theta [augmentation]
%  7  u          force applied to cart [action]
%

% 1b. Important indices

% odei  indicies for the ode solver
% augi  indicies for variables augmented to the ode variables
% dyno  indicies for the output from the dynamics model and indicies to loss
% angi  indicies for variables treated as angles (using sin/cos representation)
% dyni  indicies for inputs to the dynamics model
% poli  indicies for the inputs to the policy
% difi  indicies for training targets that are differences (rather than values)

% The indicies used here correspond to the list in 1.a

odei = [1 2 3 4];            % varibles for the ode solver 
% The ODE-solver requires to know what parts of the state are used for the 
% forward dynamics. These indices are captured by odei

augi = [];                   % variables to be augmented

dyno = [1 2 3 4];            % variables to be predicted (and known to loss)
% The predictive dimensions of the dynamics GP model are stored in dyno

angi = [4];                  % angle variables
% The indices in angi indicate which variables are angles. We represent
% these angle variables in the complex plane to exploit wrapping.

dyni = [1 2 3 5 6];          % variables that serve as inputs to the dynamics GP
% Describes which variables from the auxiliary state vector in Equation (4.1) 
% are used as the training inputs for the GP dynamics model. Note that we use
% the complex representation [sin ??, cos ??] instead of ??.

poli = [1 2 3 5 6];          % variables that serve as inputs to the policy
% Describe which state variables from the auxiliary state vector are used 
% as inputs to the policy.

difi = [1 2 3 4];            % variables that are learned via differences
% difi indices are a subset of dyno and contain the indices of the state 
% variables for which the GP training targets are differences


%% 2. Set up the scenario
dt = 0.05;                         % [s] sampling time
T = 2.5;                           % [s] initial prediction horizon time
H = ceil(T/dt);                    % prediction steps (optimization horizon)
mu0 = [0 0 0 pi+0.2]';                  % initial state mean: encodes that the cart is in the middle of the track (with zero velocity) with the pendulum hanging down (with zero angular velocity).
S0 = diag(ones(1,4)*0.01.^2);   % initial state covariance
K = 1;                             % number of initial states for which the policy is learned (how many optmization of the policy is coducted based on test scenarios)
nc = 10;                           % number of RBFs for controller parametrization

%% 3. Plant structure

% This is a structure used by the toolbox to describe the pland

plant.dynamics = @cp_dynamics;                 % function handle to the function that implements the ODE for simulating the system.
plant.noise = diag(ones(1,4)*0.01.^2);            % contains the measurement noise covariance matrix. We assume the noise is zero-mean Gaussian.
plant.dt = dt;                                    % sampling time
plant.ctrl = @zoh;                                % controler is zero order hold (Other options are @foh (first-order-hold) and @lag (first-order-lag))
plant.odei = odei;                                %  indices defined earlier into the plant structure
plant.augi = augi;
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.prop = @propagated;                         % requires a function handle to the function that computes p(xt+1) from p(xt), i.e., a one-step (probabilistic) prediction. 
% In this software package, a fairly generic function called @propagated is implemented, which computes the predictive state distribution p(xt+1) and the partial derivatives 
%that are required for gradient-based policy search.

%% 4. Policy structure
policy.fcn = @(policy,m,s)conCat(@congp,@gSat,policy,m,s);  % controller structure and parametrization (defined as a function handle)
% The policy function handle is stored in policy.fcn. The nonlinear policy here is composed from (via concatenation) a 
% deterministic GP (i.e., an RBF network) and squashing function:
%   - the RBF controller (@congp, see ctrl/congp.m), which is parametrized as the mean of a GP with a squared exponential covariance function 
%   - squashing function ?? (@gSat, see <pilco root>/util/gSat.m), defined
%     as ?? (x)=u_max(9sin(x)+sin(3x))/8 gives normalization to [???u_max, ???u_max]
policy.maxU = 40;                                         % defines input constraints, limits in the squashing function


[mm ss cc] = gTrig(mu0, S0, plant.angi);                  % augment the original state by [sin ??, cos ??] by means of gTrig.m, where the indices of the angles ?? are stored in plant.angi.
% We compute a Gaussian approximation to the joint distribution p(x,sin??,cos??).

mm = [mu0; mm];                                           % append the mean, covariance, etc. with the moments of the join distribution for the augment state
cc = S0*cc; 
ss = [S0 cc; cc' ss];        

%The following lines are used to initialize the policy parameters (definition of an inital controller):

policy.p.inputs = gaussian(mm(poli), ss(poli,poli), nc)'; % These values play the role of the training inputs of a GP and correspond to the centers ci of the RBF bassis functions.
% We sample the initial locations of the centers from the initial state distribution p(x0)

policy.p.targets = 0.1*randn(nc, length(policy.maxU));    % init. policy targets (close to zero)
% These values play the role of GP training targets and are initialized to values close to zero.
policy.p.hyp = log([1.15 3.52 3.96 0.9 0.9 0.1 0.01])';              % initialize policy hyper-parameters
% These values play the role of the GP log-hyper-parameters: log-length-scales (first 5 parameters -> augmented states), log-signal-standard-deviation (signal varaince of the controller, 6th parameter), and log-noise-standard-deviation (noise variance, relative to signal variance). 
% Expalantion:
% - log-length-scales: These values largely depend on the scale of the input data. In our example, the cart position and velocity are measured in m and m/s, respectively, the angular velocity is measured in rad/s. The last two length-scales scale trigonometric values, i.e., sin(??) and cos(??). Since these trigonometric functions map their argument into the interval [???1, 1], we choose a length-scale of 0.7, which is somewhat smaller than unity.
% - log-signal-standard-deviation: We set the signal variance of the controller to 1. Note that we squash the inital policy through ??. To exploit the full potential of the squashing function ??, it is sufficient to cover the domain [?????/2, ??/2]. Therefore, we initialize the signal variance to 1.
% - The noise variance is only important only important as a relative factor to the signal variance. This ratio essentially determines how smooth the policy is. We initialize the noise variance to 0.012.

%% 5. Set up the cost structure
cost.fcn = @cp_loss;                       % cost function handle: saturating cost function (an unnormalized Gaussian subtracted from 1) [see Sol_3_loss.m]
cost.gamma = 1;                            % discount factor (as in standard RL)

%The following parameters are specific to the cost function: [see Sol_3_loss.m]
cost.p = 0.6;                              % length of pendulum, to compute the Euclidean distance of the tip of the pendulum from the desired location in the inverted position
cost.width = 0.25;                         % cost function width, rule of thumb: cost.width= (mu0 ??? xtarget)/10
cost.expl =  0.0;                          % exploration parameter (UCB), Negative values encourage explo- ration, positive values encourage the policy staying in regions with good predictive performance. We set the value to 0 in order to disable any kind of additional exploration or exploitation.
cost.angle = plant.angi;                   % indices of angles, tells the cost function, which indices of the state are angles. 
cost.target = [0 0 0 pi]';                 % target state [set-point control]

%% 6. Dynamics model structure
% Defines the inference and training structure. We generally assume that the model uses 
% a squared exponential covariance, a Gaussian likelihood, and a zero prior 
% mean. Therefore, these parameters are not explicitly specified here.

dynmodel.fcn = @gp1d;                % built in function for GP predictions, which can predict with (sparse) GPs at uncertain inputs. If the GP is not sparse but full, gp1d calls gp0d, which implements GP predictions at uncertain inputs with the full GP model.
dynmodel.train = @train;             % built in function to train dynamics model, responsible for GP training.
dynmodel.induce = zeros(300,0,1);    % is optional and tells when to switch from full GPs to sparse GPs [for advanced users]
trainOpt = [300 500];                % contains the number of line searches for GP hyper-parameter training used by the full GP (first number) and the sparse GP (second parameter).

%% 7. Parameters for policy optimization

% In the following lines, we define (optional) parameters for policy learning. 
% Generally, an adapted version of minimize.m3, a non-convex gradient-based 
% optimizer, is used. Optional parameters are stored in a structure opt:

opt.length = 150;                        % sets the maximum number of line searches after which the optimizer returns the best parameter set so far.
opt.MFEPLS = 30;                         % is the maximum number of function evaluations per line search. Either the line search succeeds by finding a parameter set with a gradient close to 0 or it does not succeed and aborts after opt.MFEPLS many function (and gradient) evaluations.
opt.verbosity = flag.printing;           % regulates the verbosity of the optimization procedure. Verbosity ranges from 0 (no information displayed) to 3 (visualize the line search and the computed gradients).

%% 8. Plotting verbosity
plotting.verbosity = flag.plotting;      % how much information is visualized during policy learning: 0: no plots; 1: some plots; 2: all plots

%% 9. Some initializations

fantasy.mean = cell(1,N_episode);                % Used in learnPolicy to store predicted cost trajectory. This way we can keep track how predictions matched with realzied cost
fantasy.std = cell(1,N_episode);
realCost = cell(1,N_episode);                    % Used to store achieved cost trajectory during interaction with the system. Results of each episode are stored here
M = cell(N_episode,1);                           % Used to store the mean and variance of the predicted state tube per episode
Sigma = cell(N_episode,1);


mu0Sim(odei,:) = mu0;       % Define that during simulation (prediction in policy optimization) the inital state distribution is the same as during experiments
S0Sim(odei,odei) = S0;
mu0Sim = mu0Sim(dyno);      % reduce it to variables to be predicted
S0Sim = S0Sim(dyno,dyno);

