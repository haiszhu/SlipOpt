% check nullspace
% 
% 1. check rank for SLP operator, SLP traction operator, and operator that
%    combines the tangential of traction and normal of single layer.
%    expecting ~50% rank deficiency, not sure about the mixed BVP
%
% 2. check normal density contribution to above three operators
%    expecting 0 contribution
%
% 3. check baseline GMRES error
%    expecting comparable performance with preconditioner
%
% 4. check SVD pseduo-inverse solve error
%    what we are currently use, simplest setup
%
% 09/22/26 Hai

close all
clear all

root_dir = pwd;
if ~isfolder(fullfile(root_dir, 'matlab'))
  root_dir = fileparts(root_dir);
end
addpath(fullfile(root_dir, 'matlab'));
addpath(fullfile(root_dir, 'src'));
addpath(fullfile(root_dir, 'test'));

p_list = [12, 16, 20, 24, 28, 32];
tol_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
% tol_list = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4];
% p_list = [24];
% tol_list = [1e-5];
results = nan(numel(p_list), 8);

for ip = 1:numel(p_list)

  p = p_list(ip);

  % p = 16;
  % p = 24;
  % % p = 32;

  Shape = 'Y43';
  factor = 1.5;
  
  nu = p + 1;
  nv = 2 * p;
  np = nu * nv;
  r = zeros(3, np);
  nx = zeros(3, np);
  xu = zeros(3, np);
  xv = zeros(3, np);
  w = zeros(np, 1);
  [r, nx, xu, xv, w] = build_geometry_bvp_mex(p, np, Shape, r, nx, xu, xv, w);
  fprintf('Building matrices for shape=%s, p=%d (size %d x %d)\n', Shape, p, 3*np, 3*np);
  
  %% build SMat and TMat, and check rank
  % the rank should be aligned with num of spherical-harmonics vs nodal representation
  SMat = zeros(3*np, 3*np);
  SMat = stokesKernelMat_mex(p, np, 3*np, r, w(:), SMat);
  TMat = zeros(3*np, 3*np);
  TMat = kerneldSMatrix_mex(p, np, 3*np, r, nx, w(:), TMat);
  TractionMat = TMat - 0.5*eye(3*np);
  MixMat = zeros(3*np);
  t1 = xu ./ vecnorm(xu);
  t2 = xv ./ vecnorm(xv);
  for c = 1:3
    MixMat(1:3:end,:) = MixMat(1:3:end,:) + t1(c,:)' .* TractionMat(c:3:end,:);
    MixMat(2:3:end,:) = MixMat(2:3:end,:) + t2(c,:)' .* TractionMat(c:3:end,:);
    MixMat(3:3:end,:) = MixMat(3:3:end,:) + nx(c,:)' .* SMat(c:3:end,:);
  end
  % spherical-harmonic modes: \sum_0^p (2*ell+1) = (p+1)^2; 
  % among all, one of sin(ell*2*pi/(2*p)) = 0, which gives rank (p+1)^2-1
  sh_rank_bound = 3*((p+1)^2 - 1);
  [US, DS, VS] = svd(SMat, 'econ');
  [UK, DK, VK] = svd(TMat, 'econ');
  [UT, DT, VT] = svd(TractionMat, 'econ');
  [UM, DM, VM] = svd(MixMat, 'econ');
  sS = diag(DS);
  sK = diag(DK);
  sT = diag(DT);
  sM = diag(DM);
  tolS = max(size(SMat)) * eps(sS(1));
  tolK = max(size(TMat)) * eps(sK(1));
  tolT = max(size(TMat)) * eps(sT(1));
  tolM = max(size(MixMat)) * eps(sM(1));
  rankS = nnz(sS > tolS);
  rankK = nnz(sK > tolK);
  rankT = nnz(sT > tolT);
  rankM = nnz(sM > tolM);
  fprintf('\nshape=%s, p=%d, N=%d\n', Shape, p, np);
  fprintf('Vector spherical-harmonic rank bound: %d\n\n', sh_rank_bound);
  fprintf('%-24s %8s %8s %8s %16s\n', 'Operator', 'Size', 'Rank', 'Nullity', 'sigma(rank)');
  fprintf('%-24s %8d %8d %8d %16.6e\n', 'SLP: S', 3*np, rankS, 3*np-rankS, sS(rankS));
  fprintf('%-24s %8d %8d %8d %16.6e\n', 'Integral kernel: K', 3*np, rankK, 3*np-rankK, sK(rankK));
  fprintf('%-24s %8d %8d %8d %16.6e\n', 'Traction: -I/2 + K', 3*np, rankT, 3*np-rankT, sT(rankT));
  fprintf('%-24s %8d %8d %8d %16.6e\n', 'Mixed operator', 3*np, rankM, 3*np-rankM, sM(rankM));
  
  %% property of Stokes SLP operator: S
  % Pozrikidis, p79, problems 3.1.1; apply divergence thm.
  resS = norm(SMat * nx(:)) / norm(nx(:));
  fprintf('SLP normal-density residual: %.6e\n', resS);
  
  %% Normal-density test for the traction operator: -I/2 + K
  % Pozrikidis, p81, Eq. (3.2.7), K[n] = n/2.
  resT = norm((TractionMat) * nx(:)) / norm(nx(:));
  fprintf('Traction normal-density residual: %.6e\n', resT);
  
  %% Normal-density test for the mixed operator: [t1^T*(-I/2 + K); t2^T*(-I/2 + K); n^T*S]
  % follow S[n] = 0 and (-I/2 + K)[n] = 0.
  resMix = norm(MixMat * nx(:)) / norm(nx(:));
  fprintf('Mixed normal-density residual: %.6e\n', resMix);
  
  %% Mixed BVP gmres: exact boundary data from an interior point force
  x_force = [0.025; 0.05; -0.075];
  F = [1; 1/2; 1/3];
  Upt = Sto3dSLPmat_mex(np, 1, 3*np, 3, r, x_force, 1, 0, zeros(3*np, 3));
  Fpt = Sto3dSLPnmat_mex(np, 1, 3*np, 3, r, nx, x_force, 1, 0, zeros(3*np, 3));
  u_exact = reshape(Upt * F, np, 3)';
  f_exact = reshape(Fpt * F, np, 3)';
  rhsmix = [sum(t1 .* f_exact, 1);
            sum(t2 .* f_exact, 1);
            sum(nx .* u_exact, 1)];
  b = rhsmix(:);
  % Baseline GMRES with 
  % 1. rank one correction (does not seem to matter much?)
  %     Pozrikidis, p123, Eq. (4.6.12), rank one normal correction. 
  %     Also randomized method by Sifuentes, Gimbutas, and Greengard, https://arxiv.org/abs/1401.3068
  % 2. preconditioner L, converts local (t1, t2, n) residual components to Cartesian (x, y, z) components.
  gmres_tol = tol_list(ip);
  maxit = min( 1000, numel(b));
  [~, gwt] = gauss(p+1);
  wsurf = w(:) .* (pi/p * repmat(gwt(:), 2*p, 1) ./ sin(gl_grid(p)));
  q = reshape(nx .* wsurf.', [], 1) / sum(wsurf);
  z = zeros(3*np, 1);
  z(3:3:end) = 1;
  Areg = @(v) MixMat*v + z*(q'*v);
  % Areg = @(v) MixMat*v; % without rank one correction 
  % define preconditioner
  L = spalloc(3*np, 3*np, 9*np);
  for i = 1:np
    idx = 3*i-2:3*i;
    L(idx,idx) = [t1(:,i), t2(:,i), nx(:,i)]' \ eye(3);
  end
  % L = speye(3*np); % without preconditioner
  [mu_reg, flag_reg, relres_reg, iter_reg, resvec_reg] = gmres(Areg, b, [], gmres_tol, maxit, @(v) L*v);
  fprintf('GMRES: flag=%d, iterations=%d, residual=%.6e\n', flag_reg, numel(resvec_reg)-1, relres_reg);
  % evaluate velocity
  [s, t, vis] = prepare_bvp_slice(p, r, w, mu_reg, factor);
  u_num = reshape(Sto3dSLPmat_mex(t.n, s.n, 3*t.n, 3*s.n, t.x, s.x, s.w, 0, zeros(3*t.n, 3*s.n)) * s.mu, [], 3)';
  u_ref = reshape(Sto3dSLPmat_mex(t.n, 1, 3*t.n, 3, t.x, x_force, 1, 0, zeros(3*t.n, 3)) * F, [], 3)';
  errS = vecnorm(u_num - u_ref) / max(vecnorm(u_ref));
  fprintf('velocity error: max=%.6e, mean=%.6e\n', max(errS), mean(errS));
  plot_bvp_slice(vis, u_num, errS);
  gmres_error = max(errS);
  
  %% Mixed BVP SVD pseudo-inverse solve
  svd_tol = 1e-10;
  idx = sM > svd_tol;
  mu_svd = VM(:,idx) * ((UM(:,idx)' * b) ./ sM(idx));
  fprintf('SVD: retained rank=%d, residual=%.6e\n', nnz(idx), norm(MixMat*mu_svd-b)/norm(b));
  [s, t, vis] = prepare_bvp_slice(p, r, w, mu_svd, factor);
  u_num = reshape(Sto3dSLPmat_mex(t.n, s.n, 3*t.n, 3*s.n, t.x, s.x, s.w, 0, zeros(3*t.n, 3*s.n)) * s.mu, [], 3)';
  u_ref = reshape(Sto3dSLPmat_mex(t.n, 1, 3*t.n, 3, t.x, x_force, 1, 0, zeros(3*t.n, 3)) * F, [], 3)';
  errS = vecnorm(u_num - u_ref) / max(vecnorm(u_ref));
  fprintf('SVD velocity error: max=%.6e, mean=%.6e\n', max(errS), mean(errS));
  plot_bvp_slice(vis, u_num, errS);
  
  %
  results(ip,:) = [p, resS, resT, resMix, gmres_error, max(errS), flag_reg, numel(resvec_reg)-1];

end

disp(array2table(results, 'VariableNames', {'p', 'SLP_null', 'traction_null', 'mixed_null', 'GMRES_error', 'SVD_error', 'GMRES_flag', 'iterations'}));
figure(3); clf;
semilogy(results(:,1), results(:,5:6), '-o');
xlabel('p');
ylabel('Maximum relative velocity error');
legend('GMRES', 'SVD', 'Location', 'best');
grid on;

keyboard


function [s, t, vis] = prepare_bvp_slice(p, r, w, mu_gmres, factor)
s = struct('x', r, 'n', size(r, 2));
[~, gwt] = gauss(p+1);
wt = pi/p * repmat(gwt(:), 2*p, 1) ./ sin(gl_grid(p));
s.w = w(:) .* wt(:);
muS = reshape(mu_gmres, 3, [])';
s.mu = muS(:);
[theta, phi] = gl_grid(p);
theta = reshape(theta, p+1, 2*p);
phi = reshape(phi, p+1, 2*p);
theta = [pi*ones(1, 2*p); theta; zeros(1, 2*p)];
phi = [phi; phi(1:2, :)];
XCap = zeros((p+3)*2*p, 3);
for jj = 1:3
  XCap(:, jj) = sumBasis(shAna(r(jj, :)'), 'Ynm', true, theta(:), phi(:));
end
vis.xsurf = reshape(XCap(:, 1), p+3, 2*p);
vis.ysurf = reshape(XCap(:, 2), p+3, 2*p);
vis.zsurf = reshape(XCap(:, 3), p+3, 2*p);
vis.xsurf = [vis.xsurf, vis.xsurf(:, 1)];
vis.ysurf = [vis.ysurf, vis.ysurf(:, 1)];
vis.zsurf = [vis.zsurf, vis.zsurf(:, 1)];
vis.gx = linspace(-2, 2, 201);
vis.gz = linspace(-2, 2, 201);
[vis.XX, vis.ZZ] = meshgrid(vis.gx, vis.gz);
tx = [vis.XX(:), zeros(numel(vis.XX), 1), vis.ZZ(:)]';
FV = surf2patch(factor*vis.xsurf, factor*vis.ysurf, factor*vis.zsurf, 'triangles');
vis.OUT = ~inpolyhedron(FV, tx');
t.x = tx(:, vis.OUT);
t.n = size(t.x, 2);
end

function plot_bvp_slice(vis, u_num, errS)
uMag = vecnorm(u_num);
uSlice = nan(size(vis.XX));
uSlice(vis.OUT) = uMag;
figure(2); clf;
colormap(parula);
for k = 1:2
  subplot(1, 2, k);
  surf(vis.xsurf, vis.zsurf, vis.ysurf, 'EdgeColor', 'none', 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.8);
  axis equal; hold on; view(0, 90);
end
subplot(1, 2, 1);
pc = pcolor(vis.gx, vis.gz, uSlice);
set(pc, 'FaceColor', 'interp', 'LineStyle', 'none');
clim([0, max(max(uMag), eps)]);
colorbar;
title('|u| on y=0 slice');
subplot(1, 2, 2);
logErr = log10(max(errS, realmin));
scatter3(vis.XX(vis.OUT), vis.ZZ(vis.OUT), zeros(nnz(vis.OUT), 1), 10, logErr, 'filled');
errLimits = [min(logErr), max(logErr)];
if errLimits(1) == errLimits(2)
  errLimits = errLimits + [-0.5, 0.5];
end
clim(errLimits);
colorbar;
title('log10 relative error on slice');
end