-- SlipOpt/Domain.lean
-- Geometric setting for the Stokes flow problem.
-- Declares Omega (fluid domain), Gamma (boundary), normalField (outward normal), mu (viscosity)
-- as section variables used throughout the formalization.
--
-- Corresponds to: Bonnet et al. (2604.07310) §2.1 (forward flow problem setup)
--                 Masoud & Stone §3.1 (equations of motion setup)

import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.MeasureTheory.Measure.Hausdorff
import Mathlib.Topology.Basic

namespace SlipOpt

open MeasureTheory

/-- The ambient space R3 (named R3 since the Unicode superscript form is not a valid Lean 4 identifier) -/
abbrev R3 := EuclideanSpace ℝ (Fin 3)

section DomainSetup

-- Fluid domain Omega subset R3: open, bounded, Lipschitz boundary.
-- Bonnet et al. §2.1: "Omega_s subset R^d bounded domain"
variable (Ω : Set R3)

-- Boundary Gamma = frontier Omega.
-- Bonnet et al. §2.1: "frontier Omega_s = Gamma"
variable (Γ : Set R3)

-- Unit outward normal n : R3 → R3, defined a.e. on Gamma.
-- Bonnet et al. §2.1: "n denotes the unit normal on Gamma directed away from the fluid"
variable (normalField : R3 → R3)

-- Dynamic viscosity mu > 0.  Appears in Stokes equations as coefficient of Delta u.
variable (μ : ℝ)

-- Omega is open
variable (hΩ_open : IsOpen Ω)

-- Omega is bounded
variable (hΩ_bounded : Bornology.IsBounded Ω)

-- Gamma is the topological boundary of Omega
variable (hΓ : Γ = frontier Ω)

-- Viscosity is positive
variable (hμ : 0 < μ)

-- Normal field is unit length a.e. on Gamma
variable (hn_unit : ∀ x ∈ Γ, ‖normalField x‖ = 1)

end DomainSetup

end SlipOpt
