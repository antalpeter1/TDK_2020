%% Displaying experiments for the Cart-Pole system

% Based on a script by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.


%% High-Level Steps
% For each time step, plot the observed trajectory and the predicted
% means and covariances of the Cartesian coordinates of the tip of the
% pendulum

%% Code

% Loop over states in trajectory (= time steps)
for r = 1:size(xx,1)
  if exist('j','var') && ~isempty(M{j}) % If the experiment was done with an optimzied policy
    cp_display_exp(latent{j}(r,1), latent{j}(r,4), latent{j}(r,end), cost,  ...
      ['trial # ' num2str(j+J) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec'], M{j}(:,r), Sigma{j}(:,:,r));
  else % In case of random intial experiments used to gather data
     cp_display_exp(latent{i}(r,1), latent{i}(r,4), latent{i}(r,end), cost,  ...
      ['(random) trial # ' num2str(i) ', T=' num2str(H*dt) ' sec'], ...
      ['total experience (after this trial): ' num2str(dt*size(x,1)) ...
      ' sec'])
  end
  pause(dt);
end
