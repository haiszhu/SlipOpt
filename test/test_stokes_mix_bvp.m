clear; clc;

root_dir = pwd;
if ~isfolder(fullfile(root_dir, 'matlab'))
  root_dir = fileparts(root_dir);
end
addpath(fullfile(root_dir, 'matlab'));
addpath(fullfile(root_dir, 'src'));
addpath(fullfile(root_dir, 'test'));

p = 24;
Shape = 'Y43';
factor = 1.5;

nu = p + 1;
nv = 2 * p;
np = nu * nv;
r = zeros(np, 3);
nx = zeros(np, 3);
xu = zeros(np, 3);
xv = zeros(np, 3);
w = zeros(np, 1);
[r, nx, xu, xv, w] = build_geometry_bvp_mex(p, np, Shape, r, nx, xu, xv, w);
fprintf('Building matrices for shape=%s, p=%d (size %d x %d)\n', Shape, p, 3*np, 3*np);

SMat = zeros(3*np, 3*np);
SMat = stokesKernelMat_mex(p, np, 3*np, r.', w(:), SMat);
TMat = zeros(3*np, 3*np);
TMat = kerneldSMatrix_mex(p, np, 3*np, r.', nx.', w(:), TMat);

s = struct('x', r.', 'nx', nx.');
[~, gwt] = gauss(p+1); gwt = gwt';
wt = pi/p * repmat(gwt', 2*p, 1) ./ sin(gl_grid(p));
s.w = (w(:) .* wt(:))';
xu_nrm = sqrt(sum(xu.^2, 2));
xv_nrm = sqrt(sum(xv.^2, 2));
tang1 = xu ./ xu_nrm;
tang2 = xv ./ xv_nrm;
s.tang1 = tang1.';
s.tang2 = tang2.';

y_force.x = [0.1; 0.2; -0.3] / 4;
y_force.w = 1;
pt_force = [1; 1/2; 1/3];
Ahom = zeros(3*size(s.x, 2), 3*size(y_force.x, 2));
Ahom = Sto3dSLPmat_mex(size(s.x, 2), size(y_force.x, 2), ...
  3*size(s.x, 2), 3*size(y_force.x, 2), ...
  s.x, y_force.x, y_force.w(:), 0, Ahom);
Anhom = zeros(3*size(s.x, 2), 3*size(y_force.x, 2));
Anhom = Sto3dSLPnmat_mex(size(s.x, 2), size(y_force.x, 2), ...
  3*size(s.x, 2), 3*size(y_force.x, 2), ...
  s.x, s.nx, y_force.x, y_force.w(:), 0, Anhom);
% Dirichlet boundary condition
rhs = Ahom * pt_force;
rhs = reshape(rhs, [], 3)';
% Neumann boundary condition
rhsn = Anhom * pt_force;
rhsn = reshape(rhsn, [], 3)';
% Mixed boundary condition
rhsmix1 = rhsn(1,:) .* s.tang1(1,:) + rhsn(2,:) .* s.tang1(2,:) + rhsn(3,:) .* s.tang1(3,:);
rhsmix2 = rhsn(1,:) .* s.tang2(1,:) + rhsn(2,:) .* s.tang2(2,:) + rhsn(3,:) .* s.tang2(3,:);
rhsmix3 = rhs(1,:)  .* s.nx(1,:)    + rhs(2,:)  .* s.nx(2,:)    + rhs(3,:)  .* s.nx(3,:);
rhsmix = [rhsmix1; rhsmix2; rhsmix3];

TMat = -eye(size(TMat))/2 + TMat; % traction jump
% s.tang1 & s.tang2 (1st and 2nd unknown)
MixMat1 = TMat(1:3:end,:).*s.tang1(1,:)' + TMat(2:3:end,:).*s.tang1(2,:)' + TMat(3:3:end,:).*s.tang1(3,:)';
MixMat2 = TMat(1:3:end,:).*s.tang2(1,:)' + TMat(2:3:end,:).*s.tang2(2,:)' + TMat(3:3:end,:).*s.tang2(3,:)';
% s.nx (3rd unknown)
MixMat3 = SMat(1:3:end,:).*s.nx(1,:)'    + SMat(2:3:end,:).*s.nx(2,:)'    + SMat(3:3:end,:).*s.nx(3,:)';
MixMat = zeros(3*2*p*(p+1),3*2*p*(p+1));
MixMat(1:3:end,:) = MixMat1; 
MixMat(2:3:end,:) = MixMat2; 
MixMat(3:3:end,:) = MixMat3;
[U,S,V] = svd(MixMat);

diagS = diag(S);
idx = abs(diagS) > 1e-12;
Ut = U(:,idx); 
St = S(idx,idx);
Vt = V(:,idx);
invAs = Vt * inv(St) * Ut';
muS = reshape(invAs * rhsmix(:), 3, [])';

% Build a closed surface (caps + seam) for inpolyhedron masking
Xin = s.x';
Xin = Xin(:);
[mu, mv, m] = spharm_grid_size([], size(Xin, 1) / 3);
nv = size(Xin, 2);
Xin = reshape(Xin, [], 3 * nv);

XCap = zeros((m + 3) * (2 * m), 3);
[u, v] = gl_grid(m);
u = reshape(u, mu, mv);
v = reshape(v, mu, mv);
u = [pi * ones(1, mv); u; zeros(1, mv)];
v = [v; v(1:2, :)];

for jj = 1:3
  XCap(:, jj) = sumBasis(shAna(Xin(:, jj)), 'Ynm', true, u(:), v(:));
end

xsurf = reshape(XCap(:, 1), m + 3, 2 * m);
ysurf = reshape(XCap(:, 2), m + 3, 2 * m);
zsurf = reshape(XCap(:, 3), m + 3, 2 * m);

xsurf = [xsurf, xsurf(:, 1)];
ysurf = [ysurf, ysurf(:, 1)];
zsurf = [zsurf, zsurf(:, 1)];

% Slice targets on y=0 plane, masked outside inflated surface
xlimval = [-2, 2];
zlimval = [-2, 2];
nxg = 200;
nzg = 200;
gx = linspace(xlimval(1), xlimval(2), nxg + 1);
gz = linspace(zlimval(1), zlimval(2), nzg + 1);
[XX, ZZ] = meshgrid(gx, gz);

tx = [XX(:), zeros(numel(XX), 1), ZZ(:)]';
FV = surf2patch(factor * xsurf, factor * ysurf, factor * zsurf, 'triangles');
IN = inpolyhedron(FV, tx');
OUT = ~IN;

t = struct();
t.x = [XX(OUT), zeros(nnz(OUT), 1), ZZ(OUT)]';

% Evaluate numerical and reference fields on slice targets
fhomA = zeros(3*size(t.x, 2), 3*size(y_force.x, 2));
fhomA = Sto3dSLPmat_mex(size(t.x, 2), size(y_force.x, 2), ...
  3*size(t.x, 2), 3*size(y_force.x, 2), ...
  t.x, y_force.x, y_force.w(:), 0, fhomA);
fhom = fhomA * pt_force;
fhom = reshape(fhom, [], 3)';
fnumA = zeros(3*size(t.x, 2), 3*size(s.x, 2));
fnumA = Sto3dSLPmat_mex(size(t.x, 2), size(s.x, 2), ...
  3*size(t.x, 2), 3*size(s.x, 2), ...
  t.x, s.x, s.w(:), 0, fnumA);
fnum = fnumA * reshape(muS, [], 1);
fnum = reshape(fnum, [], 3)';

errS = sqrt(sum((fnum - fhom).^2, 1)) / max(sqrt(sum(fhom.^2, 1)));

% Visualize velocity magnitude and error on the same slice
uSlice = nan(size(XX));
uMag = sqrt(fnum(1, :).^2 + fnum(3, :).^2);
uSlice(OUT) = uMag;

figure(2); clf;
subplot(1,2,1);
plt = surf(xsurf, zsurf, ysurf);
plt.EdgeAlpha = 0.0;
plt.FaceColor = 0.9 * [1 1 1];
plt.FaceAlpha = 0.8;
axis equal; hold on; view(0, 90);

pc = pcolor(gx, gz, uSlice);
set(pc, 'FaceColor', 'interp', 'LineStyle', 'none');
colormap(parula);
colorbar;
title('|u| on y=0 slice');

subplot(1,2,2);
plt2 = surf(xsurf, zsurf, ysurf);
plt2.EdgeAlpha = 0.0;
plt2.FaceColor = 0.9 * [1 1 1];
plt2.FaceAlpha = 0.8;
axis equal; hold on; view(0, 90);

scatter3(t.x(1, :), t.x(3, :), t.x(2, :), 10, log10(errS), 'filled');
colorbar;
title('log10 relative error on slice');

fprintf('Slice error stats: max=%.3e, mean=%.3e\n', max(errS), mean(errS));

keyboard
