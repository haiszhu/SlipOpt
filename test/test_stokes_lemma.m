function results = test_stokes_lemma(p, Shape)
% Port of RBSM3D/tests/Demo_Lemma_verification.m using SlipOpt MEX operators.
% References follow the current Slip Optimization draft (Secs. 2.2 and 3.1).
% Run test_stokes_lemma, or results = test_stokes_lemma(24, 'Y43').
if nargin < 1, p = 24; end
if nargin < 2, Shape = 'Y43'; end
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root_dir, 'matlab'));
addpath(fullfile(root_dir, 'src'));
svd_tol = 1e-10;
rng_state = rng;
restore_rng = onCleanup(@() rng(rng_state));

%% Geometry and Stokes operators
[theta, phi] = gl_grid(p);
np = numel(theta);
[r, nx, xu, xv, w] = build_geometry_bvp_mex(p, np, Shape, zeros(3,np), zeros(3,np), zeros(3,np), zeros(3,np), zeros(np,1));
[~, gwt] = gauss(p+1);
wsurf = w(:) .* (pi/p * repmat(gwt(:), 2*p, 1) ./ sin(theta));
wvec = repelem(wsurf, 3);
t1 = xu ./ vecnorm(xu);
t2 = xv ./ vecnorm(xv);
fprintf('Building matrices for shape=%s, p=%d (size %d x %d)\n', Shape, p, 3*np, 3*np);
SMat = stokesKernelMat_mex(p, np, 3*np, r, w(:), zeros(3*np));
TMat = kerneldSMatrix_mex(p, np, 3*np, r, nx, w(:), zeros(3*np)) - 0.5*eye(3*np);
% nx points out of the body; the paper uses the opposite normal for traction.
% Thus TMat*mu is body-normal traction; -TMat*mu is the paper traction.
SMatinv = pinv(SMat, svd_tol);
Amix = zeros(3*np);
for c = 1:3
  Amix(1:3:end,:) = Amix(1:3:end,:) - t1(c,:)' .* TMat(c:3:end,:);
  Amix(2:3:end,:) = Amix(2:3:end,:) - t2(c,:)' .* TMat(c:3:end,:);
  Amix(3:3:end,:) = Amix(3:3:end,:) + nx(c,:)' .* SMat(c:3:end,:);
end
AmixMatinv = pinv(Amix, svd_tol);

%======================= Slip optimization algorithm ======================
% step (i).
% Solve the r problems (1)-(2) with v^D = v_ell^R (1<=ell<=r), obtain
% tractions f_ell^R. Set up the resistance matrix C in (9), stored as CRmat.
vlR = zeros(3*np, 6); % basis for rigid-body velocity v^R
for c = 1:3
  e = zeros(3, 1); e(c) = 1;
  vlR(c:3:end,c) = 1;
  rotation = cross(repmat(e, 1, np), r, 1);
  vlR(:,c+3) = rotation(:);
end
% build B matrix
CU = -vlR(:,1:3);
CO = -vlR(:,4:6);
B = [CU CO];
% build C matrix: force/torque constraint block, distinct from the paper resistance matrix C.
% Integrate force and torque from the paper traction -TMat*mu.
nfMat = -(vlR(:,1:3)' .* wvec') * TMat;
ntMat = -(vlR(:,4:6)' .* wvec') * TMat;
C = [nfMat; ntMat];
D = zeros(6); % Force and torque depend on mu, with no direct U/Omega term.
% Dirichlet mobility system: [SMat B; C D]*[mu; alpha] = [vS; zeros(6,1)].
Bmix = zeros(3*np, 6);
for c = 1:3
  Bmix(3:3:end,:) = Bmix(3:3:end,:) + nx(c,:)' .* B(c:3:end,:);
end
Cmix = C;
mulR = SMatinv * vlR;
flR = -TMat * mulR;
CRmat = vlR' * (wvec .* flR);
% From (8), alpha=[U;Omega] solves CRmat*alpha=-C*SMatinv*vS.
% motion(vS) returns the coefficients alpha[vS]; the paper field RvS is vlR*motion(vS).
motion = @(vS) -(CRmat \ (C * (SMatinv * vS)));

% step (ii).
% Solve the r mixed problems (16) with f_i^alpha given by (10) (1<=i<=r), obtain
% z_i given by (15) and tractions f_i. Set up z(x) = [z_1(x),...,z_r(x)].
falpha = -flR / CRmat;
rhsmix = zeros(3*np, 6);
for c = 1:3
  rhsmix(1:3:end,:) = rhsmix(1:3:end,:) + t1(c,:)' .* falpha(c:3:end,:);
  rhsmix(2:3:end,:) = rhsmix(2:3:end,:) + t2(c,:)' .* falpha(c:3:end,:);
end
Amat = [Amix Bmix; Cmix D]; % Amat*[Mumix; Umix] = [rhsmix; zeros(6)].
mat6by6 = Cmix * AmixMatinv * Bmix;
Umix = pinv(mat6by6) * (Cmix * AmixMatinv * rhsmix);
Mumix = AmixMatinv * (rhsmix - Bmix*Umix);
ui = SMat * Mumix;
zi = ui - vlR*Umix; % obtain z_i given by (15), set up z(x)
fi = -TMat * Mumix; % obtain tractions f_i
% Check the mixed conditions (16), tangency (15), and v_i^R = Rz_i in Lemma 4(i).
residual = Amat*[Mumix; Umix] - [rhsmix; zeros(6)];
mixed_residual = norm(residual(1:3*np,:), 'fro') / norm(rhsmix, 'fro');
balance = residual(3*np+1:end,:);
force_residual = max(vecnorm(balance(1:3,:)));
torque_residual = max(vecnorm(balance(4:6,:)));
normal_slip = zeros(np, 6);
for c = 1:3
  normal_slip = normal_slip + nx(c,:)' .* zi(c:3:end,:);
end
tangency_error = max(abs(normal_slip), [], 'all');
motion_z_error = max(abs(motion(zi) - Umix), [], 'all');
disp("================================================");
disp(" Full mixed system: " + (3*np+6) + " equations, " + (3*np+6) + " unknowns per RHS; SVD cutoff = " + svd_tol);
disp(" The mixed boundary condition residual is : " + mixed_residual);
disp(" The net force and torque should be zero; their norms are : " + force_residual + ", " + torque_residual);
disp(" The normal component of z should be zero; the maximum is : " + tangency_error);
disp(" The difference between alpha[z] and Umix (Lemma 4(i)) is : " + motion_z_error);
disp(" ");

% step (iii)
% Set up matrix A in (17), A_ij = <f_i,z_j>_Gamma.
% Compute y(x) = z(x)*A^{-1} as in (19), using y = zi/A.
Araw = fi' * (wvec .* zi);
A_symmetry = norm(Araw-Araw', 'fro') / norm(Araw, 'fro');
A = (Araw + Araw') / 2;
y = zi / A;
disp("================================================");
disp(" Matrix A should be symmetric and positive definite");
disp(" The relative asymmetry before symmetrization is : " + A_symmetry);
disp(" The smallest eigenvalue is : " + min(eig(A)));
disp(" ");

%%%%%% up to here, same as optimization notebook %%%%%
%% verify (20): each y(:,i) slip produces the corresponding unit rigid-body velocity
Ay = motion(y);
basis_error = max(abs(Ay-eye(6)), [], 'all');
disp("================================================");
disp(" With normalized slip basis y");
disp(" Use v_slip = y(:,i) as slip, then it should produce a unit rigid-body motion");
disp(" The difference between rigid-body motion and eye(6) is : " + basis_error);
disp(" ");

%% verify Proposition 5(i): prescribed rigid-body motion
% given any rigid body motion we wish... output effective vslip to generate the target rigid body motion
rng(12);
beta = rand(6, 1);
slipbeta = y * beta;
solnU = motion(slipbeta);
target_error = max(abs(solnU-beta));
disp("================================================");
disp(" Given aimed rigid body motion beta = [U,O]");
disp(" Use v_slip = y*beta as slip, then compute the rigid body motion should be beta");
disp(" The difference turns out to be : " + target_error);
disp(" ");
% Verify Proposition 5(ii): P(y*beta) = beta'*A^{-1}*beta, with P defined by (11).
ubeta = slipbeta + vlR*solnU;
fbeta = -TMat * (SMatinv * ubeta);
power_beta = ubeta' * (wvec .* fbeta);
power_beta_formula = beta' * (A \ beta);
power_formula_error = abs(power_beta-power_beta_formula) / max(abs(power_beta_formula), eps);
disp("================================================");
disp(" For v_slip = y*beta, the power loss should equal beta'*(A\beta)");
disp(" The relative difference turns out to be : " + power_formula_error);
disp(" ");

%% verify the null-slip decomposition from Proposition 5(i) and linearity
% For a random slip vS, alpha[vS - y*alpha[vS]] = 0, so the difference belongs to N.
nmax = 4;
rng(12);
coeffs = rand((nmax+1)^2, 2);
vsliprand = zeros(3, np);
E = sum(xu.^2, 1);
F = sum(xu.*xv, 1);
G = sum(xv.^2, 1);
det_metric = E.*G - F.^2;
for ell = 0:nmax
  m = -ell:ell;
  Y = Ynm(ell, [], theta, phi);
  phase = exp(1i * phi * m);
  P = [zeros(numel(theta),1), Y ./ phase, zeros(numel(theta),1)];
  c = sqrt((ell+(-ell:ell+1)).*(ell+1-(-ell:ell+1))) / 2;
  Yu = -(P(:,1:end-2).*c(1:end-1) - P(:,3:end).*c(2:end)) .* phase;
  Yv = Y .* (1i*m);
  for m = -ell:ell
    j = ell^2 + ell + m + 1;
    du = Yu(:,ell+m+1).';
    dv = Yv(:,ell+m+1).';
    gradY = xu .* ((G.*du-F.*dv)./det_metric) + xv .* ((E.*dv-F.*du)./det_metric);
    curlY = cross(gradY, nx, 1);
    vsliprand = vsliprand + real(coeffs(j,1)*gradY + coeffs(j,2)*curlY);
  end
end
vsliprand = vsliprand(:);
betarand = motion(vsliprand);
vslipopt = y * betarand;
betaopt = motion(vslipopt);
solndiff = motion(vsliprand-vslipopt);
null_motion_error = max(abs(solndiff));
same_motion_error = max(abs(betaopt-betarand));
disp("================================================");
disp(" With random v_slip, compute alpha = [U,O]");
disp(" Then use v_slip - y*alpha as slip belongs to Null, the rigid body motion should be 0");
disp(" Which turns out to be : " + null_motion_error);
disp(" ");

%% verify Proposition 5(iii): power loss inequality (also Lemma 4(iii))
% Compare powers at the same rigid-body velocity, using betarand in both cases.
urand = vsliprand + vlR*betarand;
uopt = vslipopt + vlR*betaopt;
frand = -TMat * (SMatinv * urand);
fopt = -TMat * (SMatinv * uopt);
power_random = urand' * (wvec .* frand);
power_optimal = uopt' * (wvec .* fopt);
power_optimal_formula = betarand' * (A \ betarand);

disp("================================================");
disp(" With random v_slip, compute power loss");
disp(" Then use v_slip_opt as slip, the rigid body motion should be the same");
disp(" The difference in rigid-body motion is : " + same_motion_error);
disp(" But the power loss turns out to be decreasing from : " + power_random + " to " + power_optimal);
disp(" The optimal power computed via Proposition 5(ii) is : " + power_optimal_formula);
disp(" ");

% results = struct('p', p, 'shape', Shape, 'svd_tol', svd_tol, ...
%                  'mixed_residual', mixed_residual, 'force_residual', force_residual, ...
%                  'torque_residual', torque_residual, 'tangency_error', tangency_error, ...
%                  'motion_z_error', motion_z_error, 'A_symmetry', A_symmetry, ...
%                  'A_min_eigenvalue', min(eig(A)), 'basis_error', basis_error, ...
%                  'target_error', target_error, 'null_motion_error', null_motion_error, ...
%                  'same_motion_error', same_motion_error, 'power_formula_error', power_formula_error, ...
%                  'power_random', power_random, 'power_optimal', power_optimal, ...
%                  'power_optimal_formula', power_optimal_formula);
end
