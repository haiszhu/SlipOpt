# Stokes nullspace convergence

Results from [`test_stokes_nullspace.m`](test_stokes_nullspace.m) for the `Y43` surface, using an interior point force as the reference solution for the mixed boundary conditions. Here `p` is the spherical-harmonic degree, with `N = 2*p*(p+1)` surface nodes.

For the first table, GMRES used the normal rank-one correction and coordinate preconditioner `L`, with tolerances `[1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]`. The SVD solve used the original mixed matrix with an absolute singular-value cutoff of `1e-10`.

| p | SLP_null | traction_null | mixed_null | GMRES_error | SVD_error | GMRES_flag | iterations |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 0.000121440195704938 | 0.000292654002438407 | 0.000269266421679574 | 0.00327942056992365 | 0.00584089183175952 | 0 | 16 |
| 16 | 3.371256363554e-05 | 9.55150585735985e-05 | 8.40476629363724e-05 | 0.000398989559409522 | 0.0050973535992488 | 0 | 23 |
| 20 | 8.07714257633616e-06 | 2.25976589488418e-05 | 2.04727793189859e-05 | 5.48770256509344e-05 | 0.000191853039944574 | 0 | 34 |
| 24 | 2.09753180248932e-06 | 5.82900107539282e-06 | 5.27348814028434e-06 | 1.21053160275257e-06 | 7.35919816754993e-06 | 0 | 44 |
| 28 | 6.28175282321912e-07 | 1.9889260637269e-06 | 1.86133616873756e-06 | 2.42194752477644e-07 | 8.35490468499124e-07 | 0 | 82 |
| 32 | 1.87249417221334e-07 | 6.09954577561461e-07 | 5.73406192776109e-07 | 2.11333342600752e-08 | 6.17374527304779e-08 | 0 | 137 |

<table>
<tr>
<td width="24%" valign="top"><img src="figures/mixed_bvp_convergence_variable_tol.png" alt="GMRES and SVD convergence with p-dependent tolerance" width="209"></td>
<td width="76%" valign="top"><img src="figures/mixed_bvp_velocity_error_slice.png" alt="Velocity magnitude and relative error on the y=0 slice" width="659"></td>
</tr>
</table>

With the same setup and a fixed GMRES tolerance `tol_list = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]`:

| p | SLP_null | traction_null | mixed_null | GMRES_error | SVD_error | GMRES_flag | iterations |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 0.000121440195704938 | 0.000292654002438407 | 0.000269266421679574 | 0.00184305640675321 | 0.00584089183175952 | 0 | 67 |
| 16 | 3.371256363554e-05 | 9.55150585735985e-05 | 8.40476629363724e-05 | 0.000120610265460473 | 0.0050973535992488 | 0 | 43 |
| 20 | 8.07714257633616e-06 | 2.25976589488418e-05 | 2.04727793189859e-05 | 5.48770256509344e-05 | 0.000191853039944574 | 0 | 34 |
| 24 | 2.09753180248932e-06 | 5.82900107539282e-06 | 5.27348814028434e-06 | 6.14625264028823e-05 | 7.35919816754993e-06 | 0 | 32 |
| 28 | 6.28175282321912e-07 | 1.9889260637269e-06 | 1.86133616873756e-06 | 6.16673092626382e-05 | 8.35490468499124e-07 | 0 | 32 |
| 32 | 1.87249417221334e-07 | 6.09954577561461e-07 | 5.73406192776109e-07 | 6.04363231422386e-05 | 6.17374527304779e-08 | 0 | 32 |

<img src="figures/mixed_bvp_convergence_fixed_tol.png" alt="GMRES and SVD convergence with fixed GMRES tolerance of 1e-4" width="207">

# Force torque system

In [`test_stokes_lemma.m`](test_stokes_lemma.m), the density and rigid-body velocity satisfy

$$
\begin{bmatrix}
A_{\mathrm{mix}} & B_{\mathrm{mix}} \\
C_{\mathrm{mix}} & 0
\end{bmatrix}
\begin{bmatrix}
\mu \\
\alpha
\end{bmatrix}
=
\begin{bmatrix}
b \\
0
\end{bmatrix},
\qquad
\alpha =
\begin{bmatrix}
U \\
\Omega
\end{bmatrix}.
$$

A related 2D force/torque-coupled formulation is given by [Guo, Zhu, and Veerapaneni (2020), Eq. (23)](https://arxiv.org/pdf/2001.05457#page=8); grouping its first three block rows and columns yields the same `[A B; C 0]` structure.

`Amix` contains the mixed boundary conditions, `Bmix` couples the rigid-body velocity, and `Cmix = [nfMat; ntMat]` imposes zero net force and torque. The six extra equations accompany six rigid-body unknowns, giving `3*N + 6` equations and unknowns.

First compute the density-block pseudoinverse, with absolute cutoff `svd_tol = 1e-10`:

```matlab
AmixMatinv = pinv(Amix, svd_tol);
```

Using this truncated inverse, set the density from the first block row and substitute into the second:

$$
\mu=A_{\mathrm{mix}}^{\dagger}(b-B_{\mathrm{mix}}\alpha),
\qquad
(C_{\mathrm{mix}}A_{\mathrm{mix}}^{\dagger}B_{\mathrm{mix}})\alpha
=C_{\mathrm{mix}}A_{\mathrm{mix}}^{\dagger}b.
$$

Solve the resulting 6-by-6 system, then recover the density:

```matlab
mat6by6 = Cmix * AmixMatinv * Bmix;
Umix = pinv(mat6by6) * (Cmix * AmixMatinv * rhsmix);
Mumix = AmixMatinv * (rhsmix - Bmix*Umix);
```

For `Y43`, `p = 24`, and `N = 1200`. Checking the original equations after the solve gives:

| Quantity | Value |
| --- | ---: |
| Relative mixed boundary-condition residual | 1.6197e-05 |
| Net force norm | 2.2496e-10 |
| Net torque norm | 1.1881e-11 |

And various lemma verification results.