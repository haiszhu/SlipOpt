# SlipOpt

SlipOpt is a Fortran + MATLAB/MEX implementation for Stokes boundary integral workflows on smooth closed surfaces (including spherical-harmonic perturbations such as `Y43`).

## Layout

- `src/`: Fortran sources (`stokes_mod.f90`, `geometry_mod.f90`, mex wrappers, etc.)
- `matlab/`: generated and hand-written MATLAB/MEX interface files
- `test/`: end-to-end MATLAB test scripts
- `external/`: external dependencies (including `rotgrid`)

## Build

From this directory:

```bash
make
```

This builds `matlab/stokes_mex.mexmaca64` and regenerates MATLAB wrappers from `matlab/stokes.mw`.

## Core MEX APIs

- `stokesKernelMat_mex`: Stokes SLP self matrix
- `kerneldSMatrix_mex`: traction kernel self matrix
- `Sto3dSLPmat_mex`: dense Stokes SLP matrix (A output)
- `Sto3dSLPnmat_mex`: dense Stokes traction matrix (same role as `Sto3dSLPmat` second output)
- `build_geometry_bvp_mex`: geometry builder with shape name input (`'sphere'`, `'Y43'`)
- `sht_ana_rotgrid_mex`, `ynm_all_mex`: harmonic helper routines

## Geometry

Use:

```matlab
np = (p+1)*(2*p);
r  = zeros(np,3); nx = zeros(np,3);
xu = zeros(np,3); xv = zeros(np,3); w = zeros(np,1);
[r, nx, xu, xv, w] = build_geometry_bvp_mex(p, np, Shape, r, nx, xu, xv, w);
```

`build_geometry_bvp_mex` is the geometry entry point.

## Tests

Run from the `test/` folder (or project root):

- `test_stokes_dirichlet_bvp.m`
- `test_stokes_neumann_bvp.m`
- `test_stokes_mix_bvp.m`

These scripts:

- build geometry and system matrices,
- solve with truncated SVD,
- evaluate on slice targets,
- visualize magnitude/error fields.

## References

This code supports:

- Das, Kausik; Zhu, Hai; Bonnet, Marc; Veerapaneni, Shravan.
  *Squirmers with arbitrary shape and slip: modeling, simulation, and optimization*.
  [arXiv:2602.19336](https://arxiv.org/abs/2602.19336)
- Bonnet, Marc; Das, Kausik; Veerapaneni, Shravan; Zhu, Hai.
  *Slip optimization on arbitrary 3D microswimmers: a reduced-dimension and boundary-integral framework*.
  [arXiv:2604.07310](https://arxiv.org/abs/2604.07310)

This code depends on:

- `rotgrid` (Fortran)
- MATLAB files from [`ves3d/matlab`](https://github.com/mmorse1217/ves3d/tree/master/matlab) for computing the Stokes single-layer self-evaluation matrix

Parts of this codebase were co-developed with AI coding assistants, including Codex and Claude Code.
