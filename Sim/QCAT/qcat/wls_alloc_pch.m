function [u, W, iter] = wls_alloc(B, v, umin, umax, Wv, Wu, ud, gam, u, W, imax, tol, phi)

% WLS_ALLOC - Control allocation using weighted least squares with Pseudo Control Hedging.
%
%  [u, W, iter] = wls_alloc(B, v, umin, umax, [Wv, Wu, ud, gamma, u0, W0, imax])
%
% Solves the weighted, bounded least-squares problem
%
%   min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
%
%   subj. to  umin <= u <= umax
%
% using an active set method and incorporating Pseudo Control Hedging.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m) [I]
% ud    desired control (m x 1) [0]
% gamma weight (scalar) [1e6]
% u0    initial point (m x 1)
% W0    initial working set (m x 1) [empty]
% imax  max no. of iterations [100]
%
%  Outputs:
%  -------
% u     optimal control
% W     optimal active set
% iter  no. of iterations (= no. of changes in the working set + 1)
%
%                            0 if u_i not saturated
% Working set syntax: W_i = -1 if u_i = umin_i
%                           +1 if u_i = umax_i
%
% See also: WLSC_ALLOC, IP_ALLOC, FXP_ALLOC, QP_SIM.

  % Number of variables
  m = length(umin);

% Set default values of optional arguments
if nargin < 13, phi = 0;          end
if nargin < 12, tol = 1e-6;       end
if nargin < 11, imax = 100;       end
[k, m] = size(B);
if nargin < 10, u = (umin + umax) / 2; W = zeros(m, 1); end
if nargin < 8,  gam = 1e6;        end
if nargin < 7,  ud = zeros(m, 1); end
if nargin < 6,  Wu = eye(m);      end
if nargin < 5,  Wv = eye(k);      end

  gam_sq = sqrt(gam);
  % A = [gam_sq*Wv*B ; Wu];
  % b = [gam_sq*Wv*v ; Wu*ud];

  A = [gam_sq * Wv * B; Wu; sqrt(phi) * eye(m)];
  b = [gam_sq * Wv * v; Wu * ud; sqrt(phi) * u];
  % Initial residual.
  d = b - A*u;
  % Determine indices of free variables.
  i_free = W==0;

  % Initialize iteration count and residual_norm
  iter = 0;
  residual_norm = inf;
  tolerance = 1e-6; % Choose a suitable tolerance value

  % Iterate until the residual_norm is below the tolerance or the maximum number of iterations is reached
  while iter < imax && residual_norm > tolerance
    % Compute optimal perturbation vector p.
    A_free = A(:, i_free);
    p_free = A_free \ d;
    p = zeros(m, 1);
  % Insert perturbations from p_free into free variables.
p(i_free) = p_free;

% Check for feasibility of the new point.
u_opt = u + p;
infeasible = (u_opt < umin) | (u_opt > umax);

if ~any(infeasible(i_free))
  % If the new point is feasible, check for optimality.

  % Update point and residual.
  u = u_opt;
  d = d - A_free*p_free;

  % Calculate residual norm
  residual_norm = norm(d, 2);

  % Compute Lagrangian multipliers.
  lambda = W.*(A'*d);
  % Are all lambda non-negative?
  if all(lambda >= -eps)
    % Optimum found, exit the loop.
    break;
  end

  % Optimum not found, remove one active constraint.
  [lambda_neg, i_neg] = min(lambda);
  W(i_neg) = 0;
  i_free(i_neg) = 1;

else
  % If the new point is infeasible, find primary bounding constraint.

  % Compute distances to the different boundaries.
  dist = ones(m, 1);
  i_min = i_free & p < 0;
  i_max = i_free & p > 0;

  dist(i_min) = (umin(i_min) - u(i_min)) ./ p(i_min);
  dist(i_max) = (umax(i_max) - u(i_max)) ./ p(i_max);

  % Proportion of p to travel
  [alpha, i_alpha] = min(dist);
  % Update point and residual.
  u = u + alpha * p;
  d = d - A_free * alpha * p_free;

  % Add corresponding constraint to the working set.
  W(i_alpha) = sign(p(i_alpha));
  i_free(i_alpha) = 0;

end

% Increment iteration count.
iter = iter + 1;
end
end

