# An oblique piston shock and injected test species

Status: Implemented and gated (fef8bde, 2026-11-15); not merged
Refs: PIC-D12, PIC-D18 (amended: `external_fields` together with `B0` / `E0`); PIC-D13, PIC-D9 (amended: test
species enter through an injecting face); PIC-D26 (used unchanged)

## Context

The research topic `qpar-electron-conditions` (branch `research/qpar-electron-conditions`, plan approved by the
owner on 2026-11-14) runs quasi-parallel shocks: a flow along x into a reflecting wall at x−, with a magnetic
field at 8–30° to x, plasma entering through the open face at x+, and test electrons (seed electrons) that enter
with the plasma. Step 4 of that plan checks the engine with short pilots and, where the engine falls short, names
this route: a branch from main, one pull request, the GPU gate, a Codex review with no findings, the owner's
instruction to merge.

What the pilots found (main fa09726, 1 GPU, 1500×2×2 cells, mass ratio 25; scratch runs, diagnostic only):
1. **A `B0` normal to a reflecting face breaks div B and is not refused.** `grid.apply_bc_B` zeroes the normal B on
   the wall plane, so the cell at the wall has div B = B_n / dx: max |div B| = 1.1 against the contract 1e-12.
   A `B0` normal to an open face is fine (1.6e-14).
2. **The piston with an oblique field cannot be configured.** Its three parts each have one place:
   - the normal component B_x must be an external field (PIC-D12), because the wall zeroes a normal Yee B;
   - the tangential B and the motional E = −v × B must be Yee fields (`B0`, `E0`), because a wall does not short
     an external E (PIC-D18): with everything external the wall is not a piston;
   - the open face at x+ must carry the uniform incident `B0` / `E0` (PIC-D26; a table of zero amplitude does it).
   `external.check_config` refuses `external_fields` together with `B0` / `E0`. With that one refusal bypassed
   the configuration runs: Gauss ≤ 1.2e-14 and div B ≤ 2.2e-13 over 20000 steps at 30° and at 8°, and a shock
   forms (30°: speed 0.078–0.079 c and compression 3.3 against 0.0725 c and 3.48 of the MHD jump).
   The gather (`ext_sum`), the initial half kick and `total_field` already add Yee and external fields.
3. **Test species cannot enter through an injecting face.** `injection.check_config` refuses a test species in
   `injection.species`. A test species loaded in the box is swept into the wall with the flow and is not
   replaced, so the upstream holds none after one crossing time.

## Approach

**G1a — `external_fields` together with `B0` / `E0`.** Delete the refusal in `external.check_config`. Nothing else
changes: the Yee arrays hold `B0` / `E0` and their evolution, the particles feel Yee + external, and the wall and
the absorbing layer act on the Yee part only. The rule for a tangential `E0` (PIC-D18, PIC-D26) is unchanged. A
moving window still refuses both. PIC-D12 and PIC-D18 are amended.

**G1b — refuse a `B0` component normal to a reflecting face.** In `SimConfig.__post_init__` (pure Python, so the
check is tested on the CPU): a non-zero `B0[a]` with a reflecting face on axis `a` raises and names
`external_fields.B`. Not covered: an `initial_fields` entry with a B normal to a wall (unchanged here).

**G2 — test species in `injection.species`.** A test species with beams may be listed. The default
(`species: null`) stays "every plasma species with beams", so existing runs keep their bits.
- No kernel change. `inject_face` uses the charge only as `q · w` (the deposit of the cut segment, the payment
  into the bank of a paired face, the `injected_charge` sum, the partner's weight). `_inject` hands it the
  deposit charge (`qdep`, zero for a test species) instead of the force charge; for a plasma species both are
  the same number.
- A test species has no partner and no bank: it is launched as at a `flux` face at every injecting face, also at
  a `paired` one, and it leaves through a paired face as through an open face.
- `injection.check_config`: the neutraliser is a plasma species; the opposite-charge rule and the "drifts out of
  a paired face" rule apply to plasma species only; a test species without beams is refused as before.
- `injection.check_fluctuation` (PIC-D26): the charge density and the current of the Maxwell checks are those of
  the plasma species; a test species may have a `u` entry and gets its own table.
- History: `injected_kinetic` and the energy budget sum the plasma species (as `lost` and `absorbed` do);
  `injected_charge` is unchanged (a test species adds zero); `n_injected` counts every species, like `n_absorbed`.
- Ids, the random streams and the capacity already count the beams of every species.

Out of scope, recorded in the topic's notebook: the drift of the upstream field at the PIC-D26 face with the
default absorber depth (`npml = 12`; gone with `npml ≥ 24`, a config value).

## Invariants

- Gauss ≤ 1e-10 and div B ≤ 1e-12 with the combined fields: the external B has zero discrete divergence and is not
  in the Yee arrays; `B0` tangential to the wall is not touched by `apply_bc_B`.
- No new kernel `#define`, no kernel edit.
- Injection keys stay (seed, step, global cell, slot, stream): 1 GPU and N GPUs give the same run with an
  injected test species.
- No collective moves: the sums of the diagnostics change their Python-side selection only.
- Runs without `external_fields` + `B0`/`E0` and without a listed test species keep their bits.

## Verification

All criteria below are fixed here, before the first run of each test.

CPU (`tests/cpu`):
- `test_external.py`: `check_config` accepts `external_fields` with `B0` and with `E0`; the implicit solver is
  still refused. (Replaces the test of the refusal.)
- `test_external.py`: a `B0` component normal to a reflecting face raises;
  the same component normal to an open or periodic face, and a tangential one, do not. Regression: no raise
  before G1b (checked once).
- `test_injection.py`: a listed test species with beams is accepted at a flux and at a paired face; refused
  without beams and as the neutraliser; a plasma species list with one sign plus a test species of the other sign
  is still refused at a paired face (a test species is not a carrier); `check_fluctuation` accepts a table whose
  Maxwell checks hold for the plasma when a test species with a `u` entry is listed.

GPU (`tests/gpu`, 1 GPU unless noted):
- **Mixed fields give the same orbits:** in `test_external_uniform_fields_give_the_orbits_of_grid_fields`, a third
  case with part of B external and the rest in `B0` / `E0`: orbits bit-equal to the other two (FP64 and `mixed`).
  Before G1a it raises.
- **Oblique piston:** a cold flow along −x into a reflecting wall at x−, field at 30° (external B_x, `B0` B_y,
  `E0` = −v × B0), open x+ with flux injection and the zero-amplitude PIC-D26 face. Criteria: Gauss ≤ 1e-10,
  div B ≤ 1e-12; in the region farther than c t + 2 from the wall and 2 from the face, every momentum component
  within 1e-4 u of the drift and B_y, E_z within 5e-4 of their initial values (the limits of PIC-D18's test);
  the density at the wall has risen above 1.5. Before G1a it raises at construction.
- **Injected test species:** a flow through open x faces with flux injection, run without and with a listed test
  species: the plasma histories (`kinetic`, `electric`, `injected_kinetic`, `injected_charge`, `budget`,
  `n_injected` minus the test species' count) agree to relative 1e-12; Gauss ≤ 1e-10; the number of test
  particles that entered is within 5 σ of the binomial expectation of the inflow (as in
  `test_injected_flux_is_the_inflow_of_the_drifting_maxwell_juttner_beams`); ids are unique. Before G2 it raises.
- **Decomposition (≥ 4 GPUs):** `mgpu_open_flow.py` gets a variant with a listed test species (paired x−, flux y
  and z faces): ids, `n_injected`, and the histories agree between 1 GPU and a 2×2 rank grid to 1e-10, as for the
  other variants.

Change gate: `ruff check .`, `python scripts/check_kernels.py`, `pytest tests/cpu`; targeted GPU tests during the
work; one run of `pytest tests/gpu` and `python scripts/gpu_gates.py --gpus 8` on the final head, recorded in
`bench/results/`; then the Codex review until it reports no findings.

## Results

Targeted runs during the work (1 GPU, 4 GPUs for the decomposition test; diagnostic), then the full gate.

**G1a, G1b** (f539f2e):
- Split fields: the orbits of `test_external_uniform_fields_give_the_orbits_of_grid_fields` hold the same bits
  with part of each field external and the rest in `B0` / `E0`, in FP64 and `mixed`. Raises before G1a.
- `B0` normal to a reflecting face: refused; the CPU test did not raise before G1b (checked once).
- Oblique piston, as fixed above, with one addition: **the absorbing layer of the face is 24 cells deep
  (`npml = 24`).** With the default 12 the test did not meet the momentum and E_z criteria next to the face:
  momentum 2.2e-4 u (limit 1e-4), B_y 7.4e-4, E_z 1.1e-3 (limit 5e-4). The same amounts appear with an open x−
  instead of the wall, and with the whole field in `B0` (90°, no external field) the momentum is off by 5.1e-4 u:
  the deviation belongs to the PIC-D26 face at that depth, not to G1a. The criteria were not changed. With
  `npml = 24`: momentum 1.4e-5 u, B_y 4.9e-5, E_z 5.7e-5, Gauss 7.8e-15, div B 6.8e-15, wall density 2.15 times
  the loaded one. Raises before G1a.

**G2** (c329350, tests 1d5367e):
- CPU: the two tests failed before the change (the refusal of a test species; the Ampère check counted the test
  species' current: 1.49 against 0.01 of 1.45).
- Injected test species: 6615 entered against 6639 ± 77; the plasma species entered with the same counts as
  without it (6654, 6607); `kinetic`, `electric`, `injected_kinetic`, `injected_charge`, `budget` within 6.3e-16;
  Gauss 1.7e-15; ids unique. The criterion "`n_injected` minus the test species' count" is checked on the
  per-species counts at the end of the run (`sdiag`), because the history holds the total only at its records.
  Raises before G2.
- Decomposition (`seeds`): the same ids and counts on 1 GPU and a 2×2 rank grid (10873 entered; 5834 in `inject`),
  the same bank as `inject` (0.21875), histories within 1.1e-15, Gauss 2.8e-15. All six variants pass. Raises
  before G2.

**Full gate on the head fef8bde** (2026-11-15, 8 GPUs; `pytest tests/gpu` and `scripts/gpu_gates.py --gpus 8` in
`fp64`, then the gates in `mixed` and `compact`): ALL PASSED in the three precisions.
- Records: `bench/results/2026-11-15-fef8bde19-gpu-gates.json` (32 gates), `-mixed.json` (30), `-compact.json` (28).
- GPU tests: 164 passed, 0 failed, 0 skipped (161 on 9666c49; the three new ones are the oblique piston, the
  injected test species and the `seeds` variant of the decomposition test).
- Two-stream growth error: ei +3.0 %, ep −1.7 %, eip −2.3 %, ei_implicit +3.5 %, buneman −4.7 %, ei_b0 +2.2 %
  (`mixed`: buneman −4.2 %; `compact`: buneman −4.7 %, on 9666c49 −5.1 %; the rest equal).
- Gauss max over all gates 7.6e-12 (`two-stream:buneman`; `compact` 2.1e-13, `filamentation:pair`); implicit
  energy error 1.9e-10.
- Codex review of fef8bde (round 1): no findings.
