# Circuit parameters

`l4_P{2..10}_layer_shared.npz` -- layer-shared eSWAP angles of the 4x4 open
Heisenberg circuit, `theta[block, matching]` with `matching_order = (yo, ye,
xo, xe)`, `eswap(theta) = exp(-i theta SWAP / 2)`, initial state the xe dimer
covering.  P = 2..4: SCALUQ Adam (complex64, 2026-08); P = 5..10: grown one
block at a time from the previous optimum (new block at eswap(0) = I) and
re-optimised with L-BFGS-B in complex128 (lg, 2026-09-04/10).  Eight
sigma = 0.2 rad restarts at P = 6..10 return the same optimum.  Energies
(E0 = -9.189207): P=2 -8.774551, 3 -9.068455, 4 -9.126463, 5 -9.163229,
6 -9.178168, 7 -9.184673, 8 -9.186925, 9 -9.187885, 10 -9.188364.

`l6_P3_reflection_orbit_lanczos.npz` -- the 6x6, P=3 angles after the direct
fidelity optimisation with one angle per reflection bond orbit (54 angles),
started from `l4_P3_layer_shared` kicked by sigma = 0.1 rad, Adam 100 steps on
the Lanczos-corrected mVMC stream (2026-09-08).  `theta` is the flat orbit
vector, `per_layer` the number of orbits per layer (3, 6, 3, 6, ... for yo, ye,
xo, xe), `held_out_fidelity` 0.7454 on the first 5000 reference samples;
measured afterwards 0.7505 (D=32) / 0.7570 (D=48) on all 20000.

`circuit_parameters.json` -- both, as plain lists.
