# Ground-state preparation circuits

Parameters and benchmark data for Heisenberg and half-filled Hubbard ground-state preparation on open square lattices.

## Start here: parameter CSVs

All three CSVs use the same columns and retain the source angle precision.

| File | Parameters |
| --- | --- |
| `heisenberg_layer_shared.csv` | 4×4 Heisenberg, N_b = 2–10 |
| `hubbard_layer_shared.csv` | 4×4 Hubbard at U/t = 8, N_b = 3–5; includes both fixed Heisenberg theta and optimized SW phi |
| `heisenberg_reflection_orbit.csv` | Optimized 6×6 Heisenberg, N_b = 3 |

Each row is one independent angle. `L_x,L_y` specify the optimization lattice, not every lattice to which the parameters may be transferred. `N_b` counts Heisenberg blocks. `gate=eSWAP` means theta; `gate=SW` means phi. `block` is one-based within the indicated gate type: the single SW block is block 1 and is applied after all Heisenberg blocks. `angle_rad` is in radians.

Matchings `y_1,y_0,x_1,x_0` correspond to the manuscript notation and appear in application order. `sharing=layer` applies one angle to all gates in that matching; `orbit_index` is then empty. For `sharing=reflection_orbit`, `orbit_index` is the one-based position in the original layer's angle list. These files preserve that ordering; they do not supply a bond-to-orbit map.

Use the Heisenberg angles inside `hubbard_layer_shared.csv` with its SW angles to reproduce the Hubbard results. They come from different optimization records than the standalone Heisenberg CSV and must not be interchanged.

## Contents

- `heisenberg/parameters/`: layer-shared parameters optimized on the 4×4 lattice for N_b = 2–10, and the optimized 6×6 reflection-orbit parameters. See the included README for conventions and provenance.
- `hubbard/L4_parameters/`: fixed Heisenberg parameters and optimized SW charge-fluctuation angles for the 4×4 Hubbard model at U/t = 8 and N_b = 3, 4, 5, together with original optimization records and exact-fidelity results.
- `benchmarks/`: angle-plot data, optimization histories, and 4×4 PEPS convergence data copied from the manuscript visualization directory.
- `manifest.json`: original paths and SHA256 checksums for all copied files.

## Parameter conventions

Angles are in radians. The layer order is `yo, ye, xo, xe`, corresponding to the manuscript matchings B_y1, B_y0, B_x1, B_x0. The initial singlet covering is `xe`. Legacy filenames use `P` for the number of Heisenberg blocks, N_b.

The Hubbard circuit adds one SW charge-fluctuation block. Use the fixed Heisenberg angles bundled with each Hubbard result to reproduce that result: these are not necessarily identical to the separately supplied Heisenberg parameter set. The Hubbard README records the input differences and lattice-transposition convention.

## Circuit exports

`qasm/heisenberg/` and `qasm/hubbard/` contain 48 parameter-bound OpenQASM 2.0 circuits for L = 4, 6, 8, 10. Heisenberg exports cover N_b = 2–10; Hubbard exports cover N_b = 3–5 at U/t = 8. Larger lattices reuse the 4×4 angles without reoptimization. No large-lattice Hubbard fidelity is implied by the export.

The circuits start from all-zero qubits and contain singlet preparation, Heisenberg layers, and (for Hubbard) the singly occupied embedding followed by one SW block. They contain only `h`, `x`, `z`, `s`, `sdg`, `cx`, `ry`, and `rz` from `qelib1.inc`, not matrix gates. eSWAP follows the manuscript circuit diagram. Each SW gate is expanded into eight Pauli rotations and then basis changes, CNOT chains, and R_Z. An irrelevant overall phase is omitted. These are logical circuits, without hardware routing or Clifford+T approximation.

Heisenberg qubit i denotes site i; Hubbard qubits 2i and 2i+1 denote up and down occupations, with i = yL + x. `qasm/manifest.json` records parameter sources, qubit and rotation counts, and file checksums.

Regenerate with `python3 export_qasm.py` (standard library only). Run `python3 check_qasm.py` with NumPy installed to check local eSWAP/SW matrices, intermediate-site JW signs, and exported syntax/index conventions. Exported rotation counts are checked against the manuscript formulas. These checks do not constitute a new large-system fidelity calculation.

A second audit is available as `python3 audit_circuits.py` (NumPy). It reads all 48 exported files, checks parameter order, rotation angles, register sizes, counts and hashes, executes all 4×4 Heisenberg exports against direct eSWAP evolution, and compares each Hubbard circuit's pre-embedding spin state to the saved optimization state. Full 2×2 Hubbard circuits are also compared with fermionic SW evolution, including embedding and intermediate JW signs. Results are recorded in `qasm/validation.json`. In the Hubbard exports, the initial spin is held on the even qubit before the local embedding; the final occupation convention is still up=2i, down=2i+1.

The 6×6 reflection-orbit parameters remain CSV-only: an explicit bond-to-orbit mapping is required before exporting that optimized circuit.

## Provenance

The source directories are `latex_source/hubbard_L4_parameters` and `latex_source/visualize/data`. Their files are preserved unchanged here; the original directories remain available for manuscript generation. Original records may contain paths from the machines used for the calculations.
