function plot_bvp_slice(vis, u_num, errS)
uMag = vecnorm(u_num);
uSlice = nan(size(vis.XX));
uSlice(vis.OUT) = uMag;
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