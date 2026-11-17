# SETUP

## Checking Out the Repository

```
% pwd
/rd47/scratch/ASP
% git clone git@lex-gitlab.aer.com:AQAC/common_ASP.git GitCommonASP
```

This repository will be used with other model repositories and not so much on its own. Currently, common_ASP is being utilized by the ASP and STILT-ASP repositories (the master and STILT-ASP branches, respectively), and we will refer the reader to their respective READMEs for a summary on how to set up each of those models with this common set of code. Note that with ASP and STILT-ASP, common_ASP is being cloned as a submodule and is thus cloned in a "detached head" state. If common_ASP is cloned by itself, the user will start on the master branch.

# Optional settings in ASPInputDeck.in

After the last required line of `ASPInputDeck.in` (the lights-off time), a deck may give solver and evaluation
settings as `Name = value` lines, in any order, until the end of the file. Names are not case sensitive; `!` starts a
comment as elsewhere in the deck. Settings not given keep their defaults, so decks without these lines run as before.

| Setting | Default | Meaning |
|---|---|---|
| `GDTableMode` | 1 | 1 = tabulated Gibbs-Duhem integral, 0 = full DGAUS8 integral on every call |
| `GDTableCells` | 256 | table cells per reaction (16 to 65536) |
| `AqInnerSolver` | 1 | single-reaction electrolyte solver: 1 = Newton, 0 = Jacobson mass-flux bisection |
| `AqOuterSolver` | 1 | electrolyte equilibrium: 1 = coupled Newton solve, 0 = worst-first loop |
| `AqPairRule` | -1 | solid-circuit precipitation rule: 1 = on, 0 = off, -1 = on with the coupled solve (`AqOuterSolver = 1`) only |
| `AqExitFix` | 0 | worst-first loop: 0 = early exit on the last reaction's deviation, 1 = no early exit, 2 = also check every reaction at the end |
| `AqDevCache` | 1 | worst-first loop deviation cache: 1 = on, 0 = off (same results) |
| `GasEvalLinked` | 0 | gas-phase ODE and Jacobian evaluation: 0 = flat arrays, 1 = linked lists (same results) |
| `ThermoScanUp` | 0 | ThermoTest RH scan: 0 = down 0.90 to 0.20, 1 = up 0.20 to 0.90 (use a deck with initial RH 0.20) |

Example:

```
360   !! When to turn off lights in minutes
AqOuterSolver = 0    ! worst-first loop
AqPairRule    = 1
```

Older decks may give the first four settings as positional lines instead (one value per line, in the order
`GDTableMode`, `GDTableCells`, `AqInnerSolver`, `AqOuterSolver`); positional lines must come before any `Name = value`
line. The model writes every setting it uses to `Transcript.TXT` and standard output at the start of the run.

These settings replace the environment variables `ASP_FEE_OUTER` (`AqOuterSolver`), `ASP_AQ_INNER` (`AqInnerSolver`),
`ASP_FEE_PAIRRULE`, `ASP_FEE_EXITFIX`, `ASP_FEE_DEVCACHE`, `ASP_GAS_LINKED` (`GasEvalLinked`) and `ASP_THERMO_SCAN`
(`ThermoScanUp`), which the model no longer reads. The electrolyte-solver dump and replay test is in `Tests/`
(see `Tests/README.md`).
