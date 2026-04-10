# SlipOpt Plan

## Scope
This plan tracks the next validation and verification milestones for SlipOpt.

## 1. Add Lemma/Theorem Numerical Tests From Paper

### Goal
Implement reproducible numerical tests that correspond directly to key lemmas and theorems in the paper.

### Tasks
- Identify each lemma/theorem statement that has a numerical counterpart.
- Map each statement to a concrete test quantity (error metric, convergence rate, invariance, bound, etc.).
- Add MATLAB test scripts in `test/` with fixed seeds/parameters for reproducibility.
- Add pass/fail thresholds tied to the theoretical statement.
- Document each test with a short “paper reference -> test mapping” comment.

### Deliverables
- New test scripts for each selected lemma/theorem.
- A summary table (in README or a dedicated note) listing theorem/lemma ID, test file, and tolerance.

### Success Criteria
- Tests run end-to-end and reproduce expected rates/inequalities.
- Results are stable across reruns on the same build.

## 2. Add Verification vs Axisymmetric Solver (Codeveloped With Hanliang Guo)

### Goal
Cross-validate SlipOpt outputs against the axisymmetric reference solver.

### Tasks
- Define benchmark cases shared by both solvers (geometry, boundary data, physical parameters).
- Build an adapter to export/import data in matching conventions and ordering.
- Compare key outputs: matrix actions, surface quantities, velocity/traction fields, integrated observables.
- Quantify discrepancies with relative norms and max errors.
- Add regression tests for at least one “easy” and one “challenging” case.

### Deliverables
- Verification scripts and comparison plots/tables.
- A short note documenting conventions and any required normalization/sign adjustments.

### Success Criteria
- Agreement within predefined tolerances for benchmark cases.
- Any residual mismatch is explained and tracked as an issue.

## 3. Evaluate Lean 4 for Theorem-Proving Support

### Goal
Assess feasibility and value of using Lean 4 for formalizing parts of the analysis.

### Tasks
- Select a small pilot statement from the paper (minimal but nontrivial).
- Prototype formalization in Lean 4 (definitions + statement + proof sketch or full proof if feasible).
- Record effort, blockers, and tooling overhead.
- Decide whether Lean 4 is suitable for ongoing proof support.

### Deliverables
- A Lean 4 pilot folder or separate repo link.
- A decision memo: proceed / limited use / defer.

### Success Criteria
- Clear feasibility conclusion backed by pilot results.
- If proceeding, define next formalization targets and timeline.

## Suggested Execution Order
1. Lemma/theorem numerical tests
2. Axisymmetric cross-verification
3. Lean 4 pilot

## Notes
- Keep all new tests deterministic.
- Prefer mex-backed paths for performance-critical validation runs.
- Track convention differences explicitly (indexing, normalization, sign, ordering).
