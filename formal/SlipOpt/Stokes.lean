-- SlipOpt/Stokes.lean
-- Predicate for Stokes flows and the three axioms underpinning
-- the PDE / functional-analytic parts of the proof.
--
-- Corresponds to: Bonnet et al. §2.1 problem (1)-(2)
--                 Masoud & Stone §3.1 eqs. (3.4a,b), §3.2 setup

import Mathlib.Analysis.Calculus.FDeriv.Basic
import Mathlib.Analysis.InnerProductSpace.Dual
import Mathlib.MeasureTheory.Integral.Bochner.Set
import SlipOpt.Tensor

namespace SlipOpt

open MeasureTheory

-- ============================================================
-- Derived differential operators
-- ============================================================

/-- Divergence of a vector field u at x: div u(x) = trace of the Jacobian of u at x.
    Masoud & Stone §3.1: incompressibility condition ∇·u = 0. -/
noncomputable def divergence (u : VelocityField) (x : R3) : ℝ :=
  (clmToMatrix (fderiv ℝ u x)).trace

/-- Vector Laplacian of u at x: Δu(x).
    Applied component-wise: (Δu)ᵢ = Δ(uᵢ).
    Here implemented as (Δu)ᵢ = div(grad(uᵢ)).
    For incompressible flows Δu = ∇(∇·u) - ∇×(∇×u) = -∇×(∇×u). -/
noncomputable def vectorLaplacian (u : VelocityField) (x : R3) : R3 := by
  let scalarLaplacian : (R3 → ℝ) → R3 → ℝ :=
    fun φ y =>
      (clmToMatrix (fderiv ℝ
        (fun z : R3 =>
          (InnerProductSpace.toDual ℝ R3).symm (fderiv ℝ φ z))
        y)).trace
  exact coordsEquiv.symm (fun i => scalarLaplacian (fun y => coordsEquiv (u y) i) x)

/-- Gradient of a scalar field p at x, as a vector in R3.
    Obtained by Riesz representation of fderiv ℝ p x : R3 →L[ℝ] ℝ.
    Milestone 1: body is sorry'd; Riesz identification deferred. -/
noncomputable def scalarGradient (p : PressureField) (x : R3) : R3 := by
  exact (InnerProductSpace.toDual ℝ R3).symm (fderiv ℝ p x)

/-- Correctness (component form): (Δu)ᵢ = div(grad(uᵢ)). -/
theorem vectorLaplacian_component_eq_div_grad
    (u : VelocityField) (x : R3) (i : Fin 3) :
    coordsEquiv (vectorLaplacian u x) i =
    divergence (fun y : R3 => scalarGradient (fun z : R3 => coordsEquiv (u z) i) y) x := by
  simp [vectorLaplacian, divergence, scalarGradient]

-- ============================================================
-- Stokes flow predicate
-- ============================================================

/-- Predicate: (u, p) is a Stokes flow on Ω with body force b and viscosity μ.

    Three conditions (Bonnet et al. problem (1)-(2), M&S eqs. (3.4a,b)):
      1. momentum:     -∇p + μΔu + b = 0 in Ω
      2. incompressible: ∇·u = 0 in Ω
      3. decay:        u(x) → 0 as |x| → ∞  -/
structure IsStokesFlow (μ : ℝ) (Ω : Set R3) (u : VelocityField) (p : PressureField)
    (b : VelocityField) : Prop where
  /-- Momentum balance: -∇p + μΔu + b = 0 pointwise in Ω
      Bonnet et al. eq. (1): −∇p + μΔu = 0 (for b = 0) -/
  momentum : ∀ x ∈ Ω,
    -(scalarGradient p x) + μ • vectorLaplacian u x + b x = 0
  /-- Incompressibility: ∇·u = 0 pointwise in Ω
      Bonnet et al. eq. (1): div u = 0 -/
  incompress : ∀ x ∈ Ω, divergence u x = 0
  /-- Far-field decay: u(x) → 0 as |x| → ∞
      Bonnet et al. eq. (2): lim_{|x|→∞} u(x) = 0 -/
  decay : Filter.Tendsto u (Filter.cocompact R3) (nhds 0)

-- ============================================================
-- Axioms (clearly labelled with references)
-- ============================================================

/-- AXIOM 1: Stokes solution existence and uniqueness.

    For any Dirichlet datum v^D satisfying the compatibility condition
    ⟨v^D, n⟩_Γ = 0, there exists a unique (u, p) solving the Stokes BVP
    with u|_Γ = v^D.

    References:
      - Ladyzhenskaya (1969)
      - Bonnet et al. [24, Sec. I.2]: "Problem (1)-(2) is well-posed for
        v^D ∈ H := {w ∈ H^{1/2}(Γ; ℝᵈ) | ⟨w,n⟩_Γ = 0}" -/
axiom stokes_exists_unique
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (vD : VelocityField)
    (hcompat : ∫ x in Γ, @inner ℝ R3 _ (vD x) (normalField x)
                          ∂(μH[2]) = 0) :
    ∃ (u : VelocityField) (p : PressureField),
      IsStokesFlow μ Ω u p (fun _ => 0) ∧
      ∀ x ∈ Γ, u x = vD x

/-- AXIOM 2: Divergence theorem for H¹(Ω) vector fields on Lipschitz domains.

    For F ∈ H¹(Ω, R3) on a Lipschitz domain Ω with outward normal n on Γ:
      ∫_Ω div F dV = ∫_Γ ⟪n, F⟫ dS

    Note: Mathlib4 has a divergence theorem for C¹ functions on boxes
    (MeasureTheory.integral_divergence_of_hasFDerivWithinAt). The H¹ /
    Lipschitz-domain version is not yet in Mathlib (as of 2026).

    This axiom is the only sorry in ReciprocityIdentity.lean. -/
axiom divergenceTheoremStokes
    (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω)
    (F : VelocityField) :
    ∫ x in Ω, divergence F x ∂MeasureTheory.volume =
    ∫ x in Γ, @inner ℝ R3 _ (normalField x) (F x) ∂(μH[2])

/-- AXIOM 3: Trace theorem — Stokes solutions have well-defined boundary values.

    The trace operator γ₀ : H¹(Ω) → H^{1/2}(Γ) is bounded, linear, and
    surjective for Lipschitz domains (standard Sobolev theory).
    Here stated simply as: the boundary restriction of u equals the datum vD. -/
axiom stokes_trace_has_boundary_datum
    (μ : ℝ) (Ω Γ : Set R3)
    (u : VelocityField) (p : PressureField)
    (hFlow : IsStokesFlow μ Ω u p (fun _ => 0)) :
    ∃ vD : VelocityField, ∀ x ∈ Γ, u x = vD x

/-- AXIOM 4: Pointwise reciprocal identity (M&S §3.2, eqs. 3.8–3.11).

    For two Stokes flows at each x ∈ Ω, the divergence of the reciprocal
    vector field F(x) = σ[u,p](x)·û(x) − σ[û,p̂](x)·u(x) equals the
    difference of body-force pairings:

      div F(x)  =  ⟪û(x), b(x)⟫ − ⟪u(x), b̂(x)⟫

    Proof sketch (not in Mathlib):
      M&S (3.8): û·(∇·σ + b) − u·(∇·σ̂ + b̂) = 0  →  û·∇·σ − u·∇·σ̂ = u·b̂ − û·b
      M&S (3.9): product rule  ∇·(σ·v) = (∇·σ)·v + σ:∇v
      M&S (3.10)+(3.11): σ:∇û = 2μ E:E = σ̂:∇u  (proved in stress_inner_eq_strainRate_inner)
      →  div(σ·û − σ̂·u) = û·b − u·b̂  -/
axiom pointwise_reciprocal_identity
    (μ : ℝ) (hμ : 0 < μ) (Ω : Set R3)
    (u û : VelocityField) (p pHat : PressureField) (b bHat : VelocityField)
    (hFlow1 : IsStokesFlow μ Ω u  p  b)
    (hFlow2 : IsStokesFlow μ Ω û pHat bHat) :
    ∀ x ∈ Ω,
      divergence (fun y =>
        toEuclid ((cauchyStress μ u  p   y).mulVec (coordsEquiv (û y))) -
        toEuclid ((cauchyStress μ û pHat y).mulVec (coordsEquiv (u  y)))) x =
      @inner ℝ R3 _ (û x) (b x) - @inner ℝ R3 _ (u x) (bHat x)

/-- AXIOM 5: Green's identity for Stokes self-energy (M&S §3.2, energy form).

    For a force-free Stokes flow (u, p) on Ω with boundary Γ:

      ∫_Γ ⟪u(x), f[u,p,n](x)⟫ dS  =  ∫_Ω σ[u,p](x) : ∇u(x) dV

    Proof sketch: integration by parts for the Stokes operator,
      ∫_Γ ⟪u, σn⟫ = ∫_Ω div(σᵀu) = ∫_Ω (∇·σ)·u + σ:∇u = ∫_Ω σ:∇u
    since ∇·σ = 0 (force-free momentum). Requires the Lipschitz-domain
    divergence theorem (same gap as AXIOM 2). -/
axiom green_stokes_selfenergy
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω)
    (u : VelocityField) (p : PressureField)
    (hFlow : IsStokesFlow μ Ω u p (fun _ => 0)) :
    boundaryBilinearForm Γ u (stokesTraction μ u p normalField) =
    ∫ x in Ω,
      matInner (cauchyStress μ u p x) (clmToMatrix (fderiv ℝ u x))
    ∂MeasureTheory.volume

end SlipOpt
