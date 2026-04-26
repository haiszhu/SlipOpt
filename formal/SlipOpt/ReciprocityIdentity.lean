-- SlipOpt/ReciprocityIdentity.lean
-- Lemma 1 (Lorentz reciprocity identity) from Bonnet et al. (2604.07310).
-- Proof follows Masoud & Stone §3.2 step-by-step.
--
-- SORRY INVENTORY for this file:
--   reciprocal_identity_general: sorry-free (integrability as explicit hypotheses)
--   lemma1_energy:               sorry-free (MeasurableSet Ω as explicit hypothesis)
--   lemma1_symmetry:             sorry-free (integrability as explicit hypotheses)
--   lemma1_bilinear_sym_pos:     sorry-free
--   All sorry's live in Stokes.lean (AXIOM 2, 4, 5).
--
-- Equation markers use M&S = Masoud & Stone (2019) and BA = Bonnet et al. (2604.07310).

import SlipOpt.Stokes

namespace SlipOpt

open MeasureTheory

-- ============================================================
-- General reciprocal identity  (M&S §3.2, eqs. 3.8–3.12)
-- ============================================================

/-- General Lorentz reciprocal identity  (M&S eq. 3.12).

    For two Stokes flows (u, σ, b) and (û, σ̂, bHat) on domain Ω with boundary Γ:

      ∫_Γ ⟪û, f[u]⟫ dS − ∫_Γ ⟪u, f[û]⟫ dS
        = ∫_Ω ⟪û, b⟫ dV − ∫_Ω ⟪u, bHat⟫ dV

    Proof outline (Goal phase will fill these in):

      M&S (3.8):  û·(∇·σ + b) − u·(∇·σ̂ + bHat) = 0     [both Stokes equations]
                  ⟹ û·(∇·σ) − u·(∇·σ̂) = u·bHat − û·b

      M&S (3.9):  (∇·σ)·û = ∇·(σ·û) − σ:∇û            [product rule]
                  (∇·σ̂)·u = ∇·(σ̂·u) − σ̂:∇u

      M&S (3.10): σ:∇û = −p·(∇·û) + 2μ·E[u]:∇û
                        = 2μ·E[u]:E[û]                  [incompress + symmetry]

      M&S (3.11): 2μ·E[u]:E[û] is symmetric in u↔û
                  ⟹ σ:∇û = σ̂:∇u
                  ⟹ ∇·(σ·û) − ∇·(σ̂·u) = u·bHat − û·b  [pointwise in Ω]

      M&S (3.12): Integrate over Ω; apply divergenceTheoremStokes [AXIOM 2] -/
theorem reciprocal_identity_general
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω) (hΩ : MeasurableSet Ω)
    (u û : VelocityField) (p pHat : PressureField) (b bHat : VelocityField)
    (hFlow1 : IsStokesFlow μ Ω u  p  b)
    (hFlow2 : IsStokesFlow μ Ω û pHat bHat)
    -- Regularity: integrability of body-force pairings over Ω
    (hIntb    : Integrable (fun x => @inner ℝ R3 _ (û x) (b    x)) (volume.restrict Ω))
    (hIntbHat : Integrable (fun x => @inner ℝ R3 _ (u  x) (bHat x)) (volume.restrict Ω))
    -- Regularity: integrability of traction pairings over Γ
    (hBdy1 : Integrable (fun x => @inner ℝ R3 _ (û x) (stokesTraction μ u  p    normalField x))
               ((μH[2]).restrict Γ))
    (hBdy2 : Integrable (fun x => @inner ℝ R3 _ (u  x) (stokesTraction μ û pHat normalField x))
               ((μH[2]).restrict Γ)) :
    boundaryBilinearForm Γ û (stokesTraction μ u  p    normalField) -
    boundaryBilinearForm Γ u  (stokesTraction μ û pHat normalField)
    =
    (∫ x in Ω, @inner ℝ R3 _ (û x) (b    x) ∂volume) -
    (∫ x in Ω, @inner ℝ R3 _ (u  x) (bHat x) ∂volume) := by
  -- F(x) = σ[u,p](x)·û(x) − σ[û,p̂](x)·u(x)
  let F : VelocityField := fun x =>
    toEuclid ((cauchyStress μ u  p    x).mulVec (coordsEquiv (û x))) -
    toEuclid ((cauchyStress μ û pHat  x).mulVec (coordsEquiv (u  x)))
  -- M&S (3.8)–(3.11): pointwise div F = ⟪û,b⟫ − ⟪u,bHat⟫  [AXIOM 4]
  have hpt := pointwise_reciprocal_identity μ hμ Ω u û p pHat b bHat hFlow1 hFlow2
  -- Integrate: ∫_Ω div F dV = ∫_Ω ⟪û,b⟫ dV − ∫_Ω ⟪u,bHat⟫ dV
  have hint : ∫ x in Ω, divergence F x ∂volume =
      (∫ x in Ω, @inner ℝ R3 _ (û x) (b x) ∂volume) -
      (∫ x in Ω, @inner ℝ R3 _ (u x) (bHat x) ∂volume) := by
    have h1 : ∫ x in Ω, divergence F x ∂volume =
        ∫ x in Ω, (@inner ℝ R3 _ (û x) (b x) - @inner ℝ R3 _ (u x) (bHat x)) ∂volume :=
      MeasureTheory.integral_congr_ae
        (MeasureTheory.ae_restrict_of_forall_mem hΩ hpt)
    rw [h1, MeasureTheory.integral_sub hIntb hIntbHat]
  -- M&S (3.12): divergence theorem  ∫_Ω div F dV = ∫_Γ ⟪n, F⟫ dS  [AXIOM 2]
  have hdiv := divergenceTheoremStokes Ω Γ normalField hΓ F
  -- Rewrite ∫_Γ ⟪n, F⟫ dS as boundary bilinear forms
  --   ⟪n, σ·û⟫ = ⟪û, σ·n⟫  [inner_mulVec_symm + cauchyStress_symm]
  have hbdy : ∫ x in Γ, @inner ℝ R3 _ (normalField x) (F x) ∂(μH[2]) =
      boundaryBilinearForm Γ û (stokesTraction μ u  p    normalField) -
      boundaryBilinearForm Γ u  (stokesTraction μ û pHat normalField) := by
    have hpt2 : ∀ x, @inner ℝ R3 _ (normalField x) (F x) =
        @inner ℝ R3 _ (û x) (stokesTraction μ u  p    normalField x) -
        @inner ℝ R3 _ (u  x) (stokesTraction μ û pHat normalField x) := by
      intro x
      simp only [F, stokesTraction, tractionVec, inner_sub_right]
      have h1 :
          @inner ℝ R3 _ (normalField x)
            (toEuclid ((cauchyStress μ u p x).mulVec (coordsEquiv (û x)))) =
          @inner ℝ R3 _ (û x)
            (toEuclid ((cauchyStress μ u p x).mulVec (coordsEquiv (normalField x)))) :=
        inner_mulVec_symm _ (cauchyStress_symm μ u p x) (normalField x) (û x)
      have h2 :
          @inner ℝ R3 _ (normalField x)
            (toEuclid ((cauchyStress μ û pHat x).mulVec (coordsEquiv (u x)))) =
          @inner ℝ R3 _ (u x)
            (toEuclid ((cauchyStress μ û pHat x).mulVec (coordsEquiv (normalField x)))) :=
        inner_mulVec_symm _ (cauchyStress_symm μ û pHat x) (normalField x) (u x)
      simpa [h1, h2]
    simp only [boundaryBilinearForm]
    rw [← MeasureTheory.integral_sub hBdy1 hBdy2]
    exact MeasureTheory.integral_congr_ae (ae_of_all _ hpt2)
  -- Combine
  exact hbdy.symm.trans (hdiv.symm.trans hint)

-- ============================================================
-- Lemma 1, Part 1 — energy identity  (BA Lemma 1, first statement)
-- ============================================================

/-- BA Lemma 1 (energy identity).

    Any force-free Stokes solution satisfies:
      ⟨v^D, f[v^D]⟩_Γ  =  ∫_Ω 2μ · E[u]:E[u] dV

    The boundary power equals the bulk viscous dissipation.
    Positivity (⟨v^D, f[v^D]⟩_Γ ≥ 0) follows since μ > 0 and E:E ≥ 0.

    Proof: set û = u, b = bHat = 0 in reciprocal_identity_general;
    the LHS difference collapses, giving the energy identity. -/
theorem lemma1_energy
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω) (hΩ : MeasurableSet Ω)
    (vD u : VelocityField) (p : PressureField)
    (hFlow  : IsStokesFlow μ Ω u p (fun _ => 0))
    (hDatum : ∀ x ∈ Γ, u x = vD x) :
    -- BA Lemma 1, eq. (energy): ⟨v^D, f[v^D]⟩_Γ = ∫_Ω 2μ E[u]:E[u] dV
    boundaryBilinearForm Γ vD (stokesTraction μ u p normalField) =
    ∫ x in Ω, 2 * μ * matInner (strainRate u x) (strainRate u x) ∂volume := by
  -- Step 1: Replace vD by u in the boundary integral (agree on Γ via hDatum)
  have hvD : boundaryBilinearForm Γ vD (stokesTraction μ u p normalField) =
             boundaryBilinearForm Γ u  (stokesTraction μ u p normalField) := by
    unfold boundaryBilinearForm
    refine MeasureTheory.integral_congr_ae ?_
    exact MeasureTheory.ae_restrict_of_forall_mem
      (hΓ ▸ isClosed_frontier.measurableSet)
      (fun x hx => by simpa [hDatum x hx])
  -- Step 2: Apply AXIOM 5 (Green's identity for Stokes self-energy)
  have hgreen := green_stokes_selfenergy μ hμ Ω Γ normalField hΓ u p hFlow
  -- Step 3: σ:∇u = 2μ E:E  [stress_inner_eq_strainRate_inner + incompressibility]
  have hint : ∫ x in Ω,
      matInner (cauchyStress μ u p x) (clmToMatrix (fderiv ℝ u x)) ∂volume =
      ∫ x in Ω, 2 * μ * matInner (strainRate u x) (strainRate u x) ∂volume :=
    MeasureTheory.integral_congr_ae
      (MeasureTheory.ae_restrict_of_forall_mem hΩ
        (fun x hx => stress_inner_eq_strainRate_inner μ u u p x (hFlow.incompress x hx)))
  rw [hvD, hgreen, hint]

-- ============================================================
-- Lemma 1, Part 2 — symmetry  (BA Lemma 1, second statement)
-- ============================================================

/-- BA Lemma 1 (symmetry).

    For any two force-free Stokes solutions with data v₁^D and v₂^D:
      ⟨v₁^D, f[v₂^D]⟩_Γ  =  ⟨v₂^D, f[v₁^D]⟩_Γ

    Proof: set b = bHat = 0 in reciprocal_identity_general;
    the RHS is 0, so the two boundary integrals are equal. -/
theorem lemma1_symmetry
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω) (hΩ : MeasurableSet Ω)
    (u₁ u₂ : VelocityField) (p₁ p₂ : PressureField)
    (hFlow1 : IsStokesFlow μ Ω u₁ p₁ (fun _ => 0))
    (hFlow2 : IsStokesFlow μ Ω u₂ p₂ (fun _ => 0))
    -- Regularity: integrability of traction pairings over Γ
    (hBdy1 : Integrable (fun x => @inner ℝ R3 _ (u₁ x) (stokesTraction μ u₂ p₂ normalField x))
               ((μH[2]).restrict Γ))
    (hBdy2 : Integrable (fun x => @inner ℝ R3 _ (u₂ x) (stokesTraction μ u₁ p₁ normalField x))
               ((μH[2]).restrict Γ)) :
    -- BA Lemma 1, symmetry: ⟨v₁^D, f[v₂^D]⟩_Γ = ⟨v₂^D, f[v₁^D]⟩_Γ
    boundaryBilinearForm Γ u₁ (stokesTraction μ u₂ p₂ normalField) =
    boundaryBilinearForm Γ u₂ (stokesTraction μ u₁ p₁ normalField) := by
  -- Body forces are zero, so both volume integrability conditions are trivial
  have hInt0 : ∀ (v : VelocityField),
      Integrable (fun x => @inner ℝ R3 _ (v x) ((fun _ => (0 : R3)) x)) (volume.restrict Ω) := by
    intro v; simp only [inner_zero_right]; exact integrable_zero _ _ _
  have hrec :=
    reciprocal_identity_general μ hμ Ω Γ normalField hΓ hΩ u₂ u₁ p₂ p₁
      (fun _ => 0) (fun _ => 0) hFlow2 hFlow1
      (hInt0 u₁) (hInt0 u₂) hBdy1 hBdy2
  have hsub :
      boundaryBilinearForm Γ u₁ (stokesTraction μ u₂ p₂ normalField) -
      boundaryBilinearForm Γ u₂ (stokesTraction μ u₁ p₁ normalField) = 0 := by
    simpa using hrec
  exact (sub_eq_zero.mp hsub)

-- ============================================================
-- Lemma 1, Consequence — symmetric positive bilinear form
-- ============================================================

/-- BA Lemma 1 (consequence).

    The bilinear form (v₁^D, v₂^D) ↦ ⟨v₁^D, f[v₂^D]⟩_Γ is symmetric and positive.

    Symmetry:   from lemma1_symmetry.
    Positivity: ⟨v^D, f[v^D]⟩_Γ = ∫_Ω 2μ E:E dV ≥ 0
                from lemma1_energy and hμ : 0 < μ.

    BA uses this (via C = ⟨f^R, v^R⟩_Γ) to show the resistance matrix C
    is symmetric positive definite (eqs. 8–9). -/
theorem lemma1_bilinear_sym_pos
    (μ : ℝ) (hμ : 0 < μ) (Ω Γ : Set R3) (normalField : R3 → R3)
    (hΓ : Γ = frontier Ω) (hΩ : MeasurableSet Ω)
    (u : VelocityField) (p : PressureField)
    (hFlow : IsStokesFlow μ Ω u p (fun _ => 0)) :
    -- BA Lemma 1, consequence: ⟨v^D, f[v^D]⟩_Γ ≥ 0
    0 ≤ boundaryBilinearForm Γ u (stokesTraction μ u p normalField) := by
  have henergy :
      boundaryBilinearForm Γ u (stokesTraction μ u p normalField) =
      ∫ x in Ω, 2 * μ * matInner (strainRate u x) (strainRate u x) ∂volume :=
    lemma1_energy μ hμ Ω Γ normalField hΓ hΩ u u p hFlow (by intro x hx; rfl)
  have hnonnegIntegrand :
      ∀ x : R3, 0 ≤ 2 * μ * matInner (strainRate u x) (strainRate u x) := by
    intro x
    have hE : 0 ≤ matInner (strainRate u x) (strainRate u x) :=
      matInner_self_nonneg (strainRate u x)
    nlinarith
  have hnonnegIntegral :
      0 ≤ ∫ x in Ω, 2 * μ * matInner (strainRate u x) (strainRate u x) ∂volume := by
    simpa [Measure.restrict_apply] using
      (MeasureTheory.integral_nonneg (μ := MeasureTheory.volume.restrict Ω)
        (f := fun x : R3 => 2 * μ * matInner (strainRate u x) (strainRate u x))
        hnonnegIntegrand)
  simpa [henergy] using hnonnegIntegral

end SlipOpt
