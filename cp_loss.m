%% Cost function for the cart-pole scenario

% Based on a script by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.

%% Explanation of the Cart-Pole loss function:
% The loss is 1-\exp(-0.5*d^2*a),  where a>0 is a scaling parameter and  
% d is the difference between the actual and desired position of tip of 
% the pendulum. 
%
% The mean and the variance of the loss are computed by averaging over the
% Gaussian state distribution $p(x) = \mathcal N(m,s)$ with mean $m$ 
% and covariance matrix $s$. 
%
% Derivatives of these quantities are computed when desired. 
%
%   function [L, dLdm, dLds, S2] =  Sol_3_loss(cost, m, s)
%
%
% *Input arguments:*
%
%   cost            cost structure
%     .p            length of pendulum                              [1 x  1 ]
%     .width        array of widths of the cost (summed together)
%     .expl         (optional) exploration parameter
%     .angle        (optional) array of angle indices
%     .target       target state                                    [D x  1 ]
%   m               mean of state distribution                      [D x  1 ]
%   s               covariance matrix for the state distribution    [D x  D ]
%
% *Output arguments:*
%
%   L     expected cost                                             [1 x  1 ]
%   dLdm  derivative of expected cost wrt. state mean vector        [1 x  D ]
%   dLds  derivative of expected cost wrt. state covariance matrix  [1 x D^2]
%   S2    variance of cost                                          [1 x  1 ]
%
%

%% High-Level Steps
% 1. Argument check
% 2. Precomputations
% 3. Define static penalty as distance from target setpoint
% 4. Trigonometric augmentation    
% 5. Calculate loss and derivatives 
% 6. Normalization

function [L, dLdm, dLds, S2] = cp_loss(cost, m, s)
%% 1. Argument check
% If not all parameters are defined then set them to default values

if isfield(cost,'width'); cw = cost.width; else cw = 1; end
if ~isfield(cost,'expl') || isempty(cost.expl); b = 0; else b =  cost.expl; end

%% 2.Precomputations
D0 = size(s,2);                           % state dimension
D1 = D0 + 2*length(cost.angle);           % augmented state dimension (with sin/cos)

M = zeros(D1,1);                          % initialize the output arguments to 0
M(1:D0) = m; 
S = zeros(D1); 
S(1:D0,1:D0) = s;
Mdm = [eye(D0); zeros(D1-D0,D0)]; 
Sdm = zeros(D1*D1,D0);
Mds = zeros(D1,D0*D0); 
Sds = kron(Mdm,Mdm);

%% 3. Define static penalty as distance from target setpoint
ell = cost.p;       % pendulum length
Q = zeros(D1);      % the matrix Q is computed, such that (x − x_target)'Q(x− x_target) is the squared Euclidean distance between the tip of the pendulum in the current state and the target state.
Q([1 D0+1],[1 D0+1]) = [1 ell]'*[1 ell]; 
Q(D0+2,D0+2) = ell^2;

%% 4. Trigonometric augmentation
if D1-D0 > 0  % This block is only executed if angles are present 

  target = [cost.target(:); gTrig(cost.target(:), 0*s, cost.angle)]; % The target state is augmented
    
  % State distributions are also agumented
  i = 1:D0; k = D0+1:D1;
  [M(k) S(k,k) C mdm sdm Cdm mds sds Cds] = gTrig(M(i),S(i,i),cost.angle);
  
  % Compute derivatives (for augmentation)
  X = reshape(1:D1*D1,[D1 D1]); XT = X';              % vectorized indices
  I=0*X; I(i,i)=1; ii=X(I==1)'; I=0*X; I(k,k)=1; kk=X(I==1)';
  I=0*X; I(i,k)=1; ik=X(I==1)'; ki=XT(I==1)';

  Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:);                    % chainrule
  Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
  Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
  Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
  dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
  dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);

  S(i,k) = S(i,i)*C; S(k,i) = S(i,k)';                      % off-diagonal
  SS = kron(eye(length(k)),S(i,i)); CC = kron(C',eye(length(i)));
  Sdm(ik,:) = SS*dCdm + CC*Sdm(ii,:); Sdm(ki,:) = Sdm(ik,:);
  Sds(ik,:) = SS*dCds + CC*Sds(ii,:); Sds(ki,:) = Sds(ik,:);
end

%% 5. Calculate loss

%After all the pre-computations, the expected cost is finally computed:

L = 0; 
dLdm = zeros(1,D0); 
dLds = zeros(1,D0*D0); 
S2 = 0;

for i = 1:length(cw)                    % scale mixture of immediate costs
  cost.z = target; 
  cost.W = Q/cw(i)^2;
  
  % For all widths of the cost structure, we compute the mean and the
  % variance of the saturating cost, including the derivatives with respect
  % to the mean and the covariance of the augmented sate:
  
  [r rdM rdS s2 s2dM s2dS] = lossSat(cost, M, S);
  
  % Compute the derivatives of the expected cost and the variance of the 
  % cost with respect to the mean and covariance matrix of the state 
  % distribution:
  L = L + r; S2 = S2 + s2;
  dLdm = dLdm + rdM(:)'*Mdm + rdS(:)'*Sdm;
  dLds = dLds + rdM(:)'*Mds + rdS(:)'*Sds;
  
  % In case of UCB exploration:
  if (b~=0 || ~isempty(b)) && abs(s2)>1e-12
    L = L + b*sqrt(s2);
    dLdm = dLdm + b/sqrt(s2) * ( s2dM(:)'*Mdm + s2dS(:)'*Sdm )/2;
    dLds = dLds + b/sqrt(s2) * ( s2dM(:)'*Mds + s2dS(:)'*Sds )/2;
  end
end

% 6. Normalization
n = length(cw); 
L = L/n; 
dLdm = dLdm/n; 
dLds = dLds/n; 
S2 = S2/n;