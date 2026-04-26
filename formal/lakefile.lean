import Lake
open Lake DSL

package «SlipOpt» where
  -- name is derived from the package identifier above

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.30.0-rc2"

@[default_target]
lean_lib «SlipOpt» where
