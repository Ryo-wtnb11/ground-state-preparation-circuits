"""Export parameter-bound OpenQASM 2.0 circuits using only the Python stdlib."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORDER = ('y_1', 'y_0', 'x_1', 'x_0')


def bonds(L, matching):
    axis, parity = matching.split('_')
    for y in range(L):
        for x in range(L):
            if axis == 'x' and x % 2 == int(parity) and x + 1 < L:
                yield y * L + x, y * L + x + 1
            if axis == 'y' and y % 2 == int(parity) and y + 1 < L:
                yield y * L + x, (y + 1) * L + x


def rotation(word, angle):
    """exp(-i angle P/2), via basis changes and a CNOT parity chain."""
    ops = []
    sites = sorted(word)
    for q in sites:
        if word[q] == 'Y':
            ops.append(('sdg', (q,), None))
        if word[q] in 'XY':
            ops.append(('h', (q,), None))
    for a, b in zip(sites, sites[1:]):
        ops.append(('cx', (a, b), None))
    ops.append(('rz', (sites[-1],), angle))
    for a, b in reversed(list(zip(sites, sites[1:]))):
        ops.append(('cx', (a, b), None))
    for q in reversed(sites):
        if word[q] in 'XY':
            ops.append(('h', (q,), None))
        if word[q] == 'Y':
            ops.append(('s', (q,), None))
    return ops


def eswap(i, j, theta):
    # Gate order in latex_source/tikz/eSWAP.tikz.
    return [('cx', (j,i), None), ('s', (j,), None), ('ry', (j,), theta/2),
            ('cx', (i,j), None), ('ry', (j,), -theta/2), ('cx', (i,j), None),
            ('sdg', (j,), None), ('x', (i,), None), ('rz', (i,), -theta/2),
            ('x', (i,), None), ('cx', (j,i), None)]


def sw(i, j, phi):
    ops = []
    for spin in (0, 1):
        a, b = 2 * i + spin, 2 * j + spin
        for extra, sign1 in ((2*i+1-spin, 1), (2*j+1-spin, -1)):
            for pa, pb, sign2 in (('X', 'Y', 1), ('Y', 'X', -1)):
                word = {q: 'Z' for q in range(a+1, b)}
                word[a], word[b] = pa, pb
                if extra in word:
                    assert word[extra] == 'Z'
                    del word[extra]
                else:
                    word[extra] = 'Z'
                ops.extend(rotation(word, sign1 * sign2 * phi / 8))
    return ops


def circuit(L, rows, hubbard):
    ops = []
    scale = 2 if hubbard else 1
    for i, j in bonds(L, 'x_0'):
        a, b = scale*i, scale*j
        ops.extend([('h', (a,), None), ('x', (b,), None),
                    ('cx', (a,b), None), ('z', (a,), None)])
    for r in rows:
        if r['gate'] != 'eSWAP':
            continue
        for i, j in bonds(L, r['matching']):
            ops.extend(eswap(scale*i, scale*j, float(r['angle_rad'])))
    if hubbard:
        for i in range(L*L):
            ops.extend([('cx', (2*i,2*i+1), None), ('x', (2*i,), None)])
        for r in rows:
            if r['gate'] == 'SW':
                for i, j in bonds(L, r['matching']):
                    ops.extend(sw(i, j, float(r['angle_rad'])))
    return ops


def export():
    manifest = []
    for model in ('heisenberg', 'hubbard'):
        source = ROOT / f'{model}_layer_shared.csv'
        with source.open() as f:
            data = list(csv.DictReader(f))
        for nb in sorted({int(r['N_b']) for r in data}):
            rows = [r for r in data if int(r['N_b']) == nb]
            for L in (4,6,8,10):
                hub = model == 'hubbard'
                ops = circuit(L, rows, hub)
                n = L*L*(2 if hub else 1)
                lines = ['OPENQASM 2.0;', 'include "qelib1.inc";',
                         f'// {model}, OBC, L={L}, N_b={nb}; angles from {source.name}',
                         '// Initial input: all qubits zero. Global phase omitted.',
                         f'qreg q[{n}];']
                for gate, qs, angle in ops:
                    parameter = f'({angle:.17g})' if angle is not None else ''
                    lines.append(f'{gate}{parameter} '+','.join(f'q[{q}]' for q in qs)+';')
                path = ROOT/'qasm'/model/f'{model}_L{L}x{L}_Nb{nb}.qasm'
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('\n'.join(lines)+'\n')
                nr = sum(g in ('rz','ry') for g,_,_ in ops)
                assert nr == (3*nb+(8 if hub else 0))*2*L*(L-1)
                assert all(0 <= q < n for _,qs,_ in ops for q in qs)
                manifest.append(dict(file=str(path.relative_to(ROOT)), L=L, N_b=nb,
                                     model=model, qubits=n, rotations=nr,
                                     source=source.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (ROOT/'qasm/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Exported {len(manifest)} circuits.')


if __name__ == '__main__':
    export()
