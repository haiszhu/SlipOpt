function [r, nx, xu, xv, w] = build_geometry_bvp_mex(p, np, shape_name, r, nx, xu, xv, w)
shape_bytes = uint8(shape_name);
shape_name_c = double([shape_bytes(:).' 0]);
slen = double(numel(shape_name_c));
mex_id_ = 'build_geometry_bvp_mex(c i int[x], c i int[x], c i char[x], c i int[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[x])';
[r, nx, xu, xv, w] = stokes_mex(mex_id_, p, np, shape_name_c, slen, r, nx, xu, xv, w, 1, 1, slen, 1, np, 3, np, 3, np, 3, np, 3, np);
end
