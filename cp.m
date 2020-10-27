

%% Using the PILCO toolbox for solving the cart pole problem


clear all; 
close all;
clc;

%% 0. Settings of the learning scenarrio

N_episode=10;  % Number of episodes
N_initial=1;   % Number of initial experiements (data gathering)
flag.printing=0; % Printing optimization information: intensity 0-3
flag.plotting=3; % Plotting  information: intensity 0-3


%% 1. Defining the system


% For each scenario, we need the following set of files, which are specific to this scenario.
% - settings.m:     A file that contains scenario-specific settings and initializations   
% - loss.m:         A cost function
% - dynamics.m:     A file that implements the ODE, which governs the dynamics
% - (optional) visualization

% First the enviroment and the actual system is specified 
cp_settings;                      % Load settings for the cart-pole system





%% 2. Initializiaton

% Here we conduct initial trials without a controller to gather data in
% terms of states, actions, and resulting responses. Trajectores have
% length H (see settings).

% When iterating with a real setup, this corresponds to random experiments
% with random/preliminary controllers

x = [];                                  % used to gather intial data
y = [];

for i = 1:N_initial
   [xx, yy, realCost{i}, latent{i}] = rollout(gaussian(mu0, S0), struct('maxU',policy.maxU), H, plant, cost); % Execute interaction with the plant using the current random policy and a given inital condition (determined by a Gaussian distribution in this case)
   x = [x; xx]; y = [y; yy];       % Gathering obtained data into x and y
   if plotting.verbosity > 0       % visualization of trajectory
     if ~ishandle(1) 
         figure(1); 
     else
         set(0,'CurrentFigure',1); 
     end
     clf(1);
     cp_display_rollout;        % plot experiment (stored in xx, yy)
   end
end
J=N_initial; % Used by the toolbox


%% 3. Policy optimization

for j = 1:N_episode % Note that the loop varaiable "j" is used in the "display_rollout" script
    
  trainDynModel;   % Step 1. train GP model of the process dynamics [off line]
  % The script that takes care of training the GP executes the following high-level steps:
  %   a. Extract states and control inputs to be matched form x and y  
  %   b. Define the training inputs and targets of the GP 
  %   c. Train the GP 
  %   d. Display GP hyper-parameters, the learned noise hyper-parameters, and the signal-to-noise ratios. This information is very valuable for debugging purposes.


  learnPolicy;     % Step 2. optimization of the policy based on the learnt model [off line, simulation phase]
  %   a. Learn the policy by calling minimize. This uses the whole
  %      inference toolchain to compute the gradients required by the
  %      minimization
  %   b. (optional) Plot overall optimization progress. (line search etc.)
  %   c. Long-term prediction of a state trajectory from p(x0) using the 
  %      learned policy by calling pred. This prediction is equivalent to 
  %      the last predicted trajectory during policy learning, i.e., 
  %      the predicted state trajectory that belongs to the the learned controller.
  %   d. The predicted state trajectory is used to compute the corresponding 
  %      distribution over immediate costs by calling calcCost. => stored
  %      in "fantasy"
  %   c. (optional) Plot the predicted immediate cost distribution as a 
  %      function of the time steps. 
  
  
  applyController; % Step 3. apply controller to the system and gather data [on line, experiment phase]
  % The specific steps of this generic script are
  %     (a) determine start state 
  %     (b) generate rollout 
  %         i. compute control signal ??(xt) [this uses the policy object]
  %         ii. simulate dynamics (or apply control to real robot) [this will use the function handle that describes the dynamics]
  %         iii. transition to state xt+1 to continue the interaction
  % This script also plots the response using figures already set up by the
  % display scripts. Furthermore, it also stores the whole trajectory
  % information in a saved file.
  
  
  % When iteracting with a real setup, Step 3 means that with the control
  % policy a full experiment is run. The rest of the steps remain
  % unchanged.
  
% Plotting the most recent trajectory (stored in xx, yy)
   disp(['controlled trial # ' num2str(j)]);
   if plotting.verbosity > 0;      % visualization of trajectory
     if ~ishandle(1); figure(1); else set(0,'CurrentFigure',1); end; clf(1);
     cp_display_rollout;
   end
end
