-- SlipOpt/Tensor.lean
-- Tensor fields for Stokes flow: strain rate E[u], Cauchy stress σ[u,p],
-- traction f[u,p,n], and boundary bilinear form ⟨·,·⟩_Γ.
--
-- Corresponds to: Masoud & Stone §3.2 eqs. (3.2), (3.10), (3.11)
--                 Bonnet et al. §2.1 traction field f[v^D]

import Mathlib.Analysis.Calculus.FDeriv.Basic
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Measure.Hausdorff
import Mathlib.Analysis.Normed.Lp.MeasurableSpace
import SlipOpt.Domain

namespace SlipOpt

open MeasureTheory Matrix

-- ============================================================
-- Type aliases
-- ============================================================

/-- A velocity field: maps points in R3 to velocity vectors in R3.
    Bonnet et al. §2.1: u, v^D, v^s ∈ H -/
abbrev VelocityField := R3 → R3

/-- A pressure field: maps points in R3 to scalar pressure values.
    Appears in Cauchy stress σ = -pI + 2μE. -/
abbrev PressureField := R3 → ℝ

-- ============================================================
-- Coordinate conversion helper
-- ============================================================

/-- Equivalence between EuclideanSpace ℝ (Fin 3) and (Fin 3 → ℝ).
    Used to extract components for matrix operations. -/
noncomputable def coordsEquiv : R3 ≃ (Fin 3 → ℝ) :=
  WithLp.equiv 2 (Fin 3 → ℝ)

/-- Extract the Jacobian matrix of a continuous linear map L : R3 →L[ℝ] R3.
    J_{ij} = i-th component of L(eⱼ), where eⱼ is the j-th standard basis vector. -/
noncomputable def clmToMatrix (L : R3 →L[ℝ] R3) : Matrix (Fin 3) (Fin 3) ℝ :=
  Matrix.of (fun i j =>
    coordsEquiv (L (coordsEquiv.symm (Pi.single j 1))) i)

-- ============================================================
-- Tensor field definitions
-- ============================================================

/-- Strain rate tensor: E[u](x) = ½(∇u(x) + (∇u(x))ᵀ)
    Masoud & Stone eq. (3.2): E = ½[∇u + (∇u)ᵀ]
    Uses fderiv for the Jacobian; result is a symmetric 3×3 matrix. -/
noncomputable def strainRate (u : VelocityField) (x : R3) : Matrix (Fin 3) (Fin 3) ℝ :=
  let J := clmToMatrix (fderiv ℝ u x)
  (1/2 : ℝ) • (J + J.transpose)

/-- Cauchy stress tensor: σ[u,p](x) = -p(x)·I + 2μ·E[u](x)
    Masoud & Stone eq. (3.2): σ = -pI + 2μE
    This is the constitutive relation for a Newtonian incompressible fluid. -/
noncomputable def cauchyStress (μ : ℝ) (u : VelocityField) (p : PressureField)
    (x : R3) : Matrix (Fin 3) (Fin 3) ℝ :=
  -(p x) • (1 : Matrix (Fin 3) (Fin 3) ℝ) + (2 * μ) • strainRate u x

/-- Convert a Fin 3 → ℝ vector back to EuclideanSpace ℝ (Fin 3). -/
noncomputable def toEuclid (v : Fin 3 → ℝ) : R3 :=
  coordsEquiv.symm v

/-- Traction vector: f(x) = σ(x) ·ᵥ n(x)
    The stress tensor dotted with the unit outward normal gives the surface traction.
    Bonnet et al. §2.1: traction field f := -pn + 2μD[u]·n on Γ -/
noncomputable def tractionVec
    (σ : R3 → Matrix (Fin 3) (Fin 3) ℝ) (nField : R3 → R3) (x : R3) : R3 :=
  toEuclid ((σ x).mulVec (coordsEquiv (nField x)))

/-- The Stokes traction field: combines cauchyStress with tractionVec.
    Bonnet et al. notation: f[v^D](x) = cauchyStress(x) ·ᵥ n(x) -/
noncomputable def stokesTraction (μ : ℝ) (u : VelocityField) (p : PressureField)
    (nField : R3 → R3) (x : R3) : R3 :=
  tractionVec (cauchyStress μ u p) nField x

/-- Matrix Frobenius inner product A:B = trace(Aᵀ·B) = Σᵢⱼ Aᵢⱼ·Bᵢⱼ
    Used in M&S (3.10): σ:∇û and E[u]:E[û] -/
noncomputable def matInner (A B : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  (A.transpose * B).trace

/-- Boundary bilinear form: ⟨v₁, f[v₂]⟩_Γ = ∫_Γ ⟪v₁(x), f₂(x)⟫ dHaus²
    Bonnet et al. (before Lemma 1): ⟨·,·⟩_Γ denotes the L²(Γ) scalar product. -/
noncomputable def boundaryBilinearForm (Γ : Set R3)
    (v₁ : VelocityField) (f₂ : VelocityField) : ℝ :=
  ∫ x in Γ, @inner ℝ R3 _ (v₁ x) (f₂ x) ∂(μH[2])

-- ============================================================
-- Key algebraic lemmas (Milestone 1: sorry'd; Goal: proved for real)
-- ============================================================

/-- M&S (3.2): Strain rate is symmetric: E[u]ᵀ = E[u]
    Follows immediately from definition as symmetrization of the Jacobian. -/
theorem strainRate_symm (u : VelocityField) (x : R3) :
    (strainRate u x).transpose = strainRate u x := by
  simp [strainRate, add_comm]

/-- Frobenius norm square nonnegativity: A:A ≥ 0. -/
theorem matInner_self_nonneg (A : Matrix (Fin 3) (Fin 3) ℝ) :
    0 ≤ matInner A A := by
  unfold matInner
  rw [Matrix.trace]
  refine Finset.sum_nonneg (fun i _ => ?_)
  rw [Matrix.diag_apply, Matrix.mul_apply]
  simp only [Matrix.transpose_apply]
  refine Finset.sum_nonneg (fun j _ => ?_)
  nlinarith

/-- M&S (3.10)+(3.11): For a Newtonian incompressible fluid,
    σ:∇û = 2μ · E[u]:E[û], which is symmetric under u ↔ û.

    Derivation:
      σ:∇û = (-pI + 2μE[u]):∇û
            = -p·trace(∇û) + 2μ·E[u]:∇û    [by linearity of Frobenius product]
            = 2μ·E[u]:E[û]                  [incompressibility: trace(∇û)=0;
                                              symmetric E has zero inner product
                                              with the skew part of ∇û]
    The RHS is symmetric in u↔û since E[u]:E[û] = E[û]:E[u]. -/
theorem stress_inner_eq_strainRate_inner
    (μ : ℝ) (u û : VelocityField) (p : PressureField) (x : R3)
    (hdiv : (clmToMatrix (fderiv ℝ û x)).trace = 0) :
    matInner (cauchyStress μ u p x) (clmToMatrix (fderiv ℝ û x)) =
    2 * μ * matInner (strainRate u x) (strainRate û x) := by
  let Ju : Matrix (Fin 3) (Fin 3) ℝ := clmToMatrix (fderiv ℝ u x)
  let Jv : Matrix (Fin 3) (Fin 3) ℝ := clmToMatrix (fderiv ℝ û x)
  have hJv : Jv.trace = 0 := by simpa [Jv] using hdiv
  have htr1 : (Ju.transpose * Jv).trace = (Ju * Jv.transpose).trace := by
    calc
      (Ju.transpose * Jv).trace = ((Ju.transpose * Jv).transpose).trace := by
        exact (Matrix.trace_transpose (Ju.transpose * Jv)).symm
      _ = (Jv.transpose * Ju.transpose.transpose).trace := by rw [Matrix.transpose_mul]
      _ = (Jv.transpose * Ju).trace := by rw [Matrix.transpose_transpose]
      _ = (Ju * Jv.transpose).trace := by exact Matrix.trace_mul_comm (Jv.transpose) Ju
  have htr2 : (Ju.transpose * Jv.transpose).trace = (Ju * Jv).trace := by
    calc
      (Ju.transpose * Jv.transpose).trace = ((Ju.transpose * Jv.transpose).transpose).trace := by
        exact (Matrix.trace_transpose (Ju.transpose * Jv.transpose)).symm
      _ = (Jv.transpose.transpose * Ju.transpose.transpose).trace := by rw [Matrix.transpose_mul]
      _ = (Jv * Ju).trace := by rw [Matrix.transpose_transpose, Matrix.transpose_transpose]
      _ = (Ju * Jv).trace := by exact Matrix.trace_mul_comm Jv Ju

  change (((-(p x)) • (1 : Matrix (Fin 3) (Fin 3) ℝ) +
      (2 * μ) • ((1 / 2 : ℝ) • (Ju + Ju.transpose))).transpose * Jv).trace
    = 2 * μ * ((((1 / 2 : ℝ) • (Ju + Ju.transpose)).transpose *
      ((1 / 2 : ℝ) • (Jv + Jv.transpose))).trace)

  simp only [Matrix.transpose_add, Matrix.transpose_smul, Matrix.transpose_transpose,
    Matrix.transpose_one, Matrix.one_mul, Matrix.mul_add, Matrix.add_mul,
    Matrix.mul_smul, Matrix.smul_mul, Matrix.trace_add, Matrix.trace_smul,
    hJv]
  rw [htr1, htr2]
  ring_nf

/-- Cauchy stress is symmetric: σᵀ = σ.
    Follows because -pI is symmetric and 2μE is symmetric (strainRate_symm). -/
theorem cauchyStress_symm (μ : ℝ) (u : VelocityField) (p : PressureField) (x : R3) :
    (cauchyStress μ u p x).transpose = cauchyStress μ u p x := by
  simp [cauchyStress, Matrix.transpose_add, Matrix.transpose_smul,
        Matrix.transpose_one, strainRate_symm]

/-- For symmetric A, ⟪v, A·w⟫ = ⟪w, A·v⟫.
    Key lemma for rewriting ∫_Γ ⟪n, σ·û⟫ as ∫_Γ ⟪û, σ·n⟫ = ∫_Γ ⟪û, f[u]⟫. -/
lemma inner_mulVec_symm (A : Matrix (Fin 3) (Fin 3) ℝ) (hA : Aᵀ = A) (v w : R3) :
    @inner ℝ R3 _ v (toEuclid (A.mulVec (coordsEquiv w))) =
    @inner ℝ R3 _ w (toEuclid (A.mulVec (coordsEquiv v))) := by
  -- Move to coordinates (`coordsEquiv`) and reduce to dot-product identities.
  rw [show v = coordsEquiv.symm (coordsEquiv v) from (coordsEquiv.symm_apply_apply v).symm,
      show w = coordsEquiv.symm (coordsEquiv w) from (coordsEquiv.symm_apply_apply w).symm]
  have hInnerL :
      @inner ℝ R3 _ (coordsEquiv.symm (coordsEquiv v))
        (toEuclid (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv w))))) =
      (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv w)))) ⬝ᵥ star (coordsEquiv v) := by
    simpa [coordsEquiv, toEuclid] using
      (EuclideanSpace.inner_toLp_toLp (coordsEquiv v)
        (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv w)))))
  have hInnerR :
      @inner ℝ R3 _ (coordsEquiv.symm (coordsEquiv w))
        (toEuclid (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv v))))) =
      (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv v)))) ⬝ᵥ star (coordsEquiv w) := by
    simpa [coordsEquiv, toEuclid] using
      (EuclideanSpace.inner_toLp_toLp (coordsEquiv w)
        (A.mulVec (coordsEquiv (coordsEquiv.symm (coordsEquiv v)))))
  rw [hInnerL, hInnerR]
  simp only [star_trivial]
  -- Now prove: (A·w)⋅v = (A·v)⋅w for symmetric A.
  calc
    A.mulVec (coordsEquiv w) ⬝ᵥ coordsEquiv v
        = coordsEquiv v ⬝ᵥ A.mulVec (coordsEquiv w) := by
            simpa using (dotProduct_comm (A.mulVec (coordsEquiv w)) (coordsEquiv v))
    _ = Matrix.vecMul (coordsEquiv v) A ⬝ᵥ coordsEquiv w := by
          simpa using (Matrix.dotProduct_mulVec (coordsEquiv v) A (coordsEquiv w))
    _ = (A.transpose.mulVec (coordsEquiv v)) ⬝ᵥ coordsEquiv w := by
          rw [Matrix.mulVec_transpose]
    _ = (A.mulVec (coordsEquiv v)) ⬝ᵥ coordsEquiv w := by
          simpa [hA]
    _ = coordsEquiv w ⬝ᵥ A.mulVec (coordsEquiv v) := by
          simpa using (dotProduct_comm (A.mulVec (coordsEquiv v)) (coordsEquiv w))
    _ = A.mulVec (coordsEquiv v) ⬝ᵥ coordsEquiv w := by
          simpa using (dotProduct_comm (coordsEquiv w) (A.mulVec (coordsEquiv v)))

end SlipOpt
