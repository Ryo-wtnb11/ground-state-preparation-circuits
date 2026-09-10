# Ground-state preparation circuits

OpenQASM circuits and optimized parameters for the Heisenberg and half-filled Hubbard models on open square lattices.

## Circuits

Ready-to-use OpenQASM 2.0 files with numerical angles:

| Directory | Model | Lattice sizes | Heisenberg blocks |
| --- | --- | --- | --- |
| `qasm/heisenberg/` | Heisenberg | 4×4, 6×6, 8×8, 10×10 | N_b = 2–10 |
| `qasm/hubbard/` | Hubbard, U/t = 8 | 4×4, 6×6, 8×8, 10×10 | N_b = 3–5 |

Larger circuits reuse the 4×4 parameters. Each circuit starts from all-zero qubits; Hubbard circuits include one SW charge-fluctuation block. Gates are decomposed into standard gates from `qelib1.inc`. Hardware routing and Clifford+T synthesis are not included.

Sites follow `i = yL + x`. Heisenberg uses qubit `i`; Hubbard uses qubits `2i` and `2i+1` for up and down occupations.

## Parameters

| CSV | Contents |
| --- | --- |
| `heisenberg_layer_shared.csv` | 4×4 Heisenberg angles |
| `hubbard_layer_shared.csv` | Fixed Heisenberg angles and optimized SW angles at U/t = 8 |
| `heisenberg_reflection_orbit.csv` | Optimized 6×6 Heisenberg angles, N_b = 3; CSV only |

Angles are in radians. `L_x,L_y` identify the optimization lattice, and `block` is one-based within each gate type. Layers are applied in the order `y_1, y_0, x_1, x_0`, starting from singlets on `x_0`. For Hubbard, use both sets of angles in its CSV together. `orbit_index` follows the source layer's orbit ordering; a bond-to-orbit map is not included.

Additional results are in `benchmarks/`, `heisenberg/parameters/`, and `hubbard/L4_parameters/`.

## Export and check

```sh
python3 export_qasm.py
python3 check_qasm.py
python3 audit_circuits.py
```

Export requires only Python; checks require NumPy. Circuit counts and checksums are in `qasm/manifest.json`, and validation results are in `qasm/validation.json`.
