"""Audit exported files against independent state evolution (NumPy only)."""
import csv
import hashlib
import json
import re
import numpy as np
from check_qasm import target_sw, check as check_local
from export_qasm import ROOT, circuit

ORDER = ('y_1','y_0','x_1','x_0')


def parse(path):
    lines=path.read_text().splitlines()
    assert lines[0]=='OPENQASM 2.0;' and lines[1]=='include "qelib1.inc";'
    n=int(re.fullmatch(r'qreg q\[(\d+)\];',lines[4])[1])
    ops=[]
    for line in lines[5:]:
        m=re.fullmatch(r'(h|x|z|s|sdg|cx|ry|rz)(?:\(([-+0-9.e]+)\))? (q\[\d+\](?:,q\[\d+\])?);',line)
        assert m,line
        g,angle,args=m.groups();qs=tuple(map(int,re.findall(r'\d+',args)))
        assert len(qs)==(2 if g=='cx' else 1) and len(set(qs))==len(qs)
        assert all(0<=q<n for q in qs)
        assert (angle is not None)==(g in ('ry','rz'))
        ops.append((g,qs,float(angle) if angle is not None else None))
    return n,ops


def evolve(ops,n):
    state=np.zeros(2**n,complex);state[0]=1
    ids=np.arange(2**n)
    fixed={'h':np.array([[1,1],[1,-1]])/np.sqrt(2),'x':np.array([[0,1],[1,0]]),
           'z':np.diag([1,-1]),'s':np.diag([1,1j]),'sdg':np.diag([1,-1j])}
    pairs={q:(ids[(ids&(1<<q))==0],ids[(ids&(1<<q))==0]|(1<<q)) for q in range(n)}
    for g,qs,a in ops:
        if g=='cx':
            state=state[ids ^ (((ids>>qs[0])&1)<<qs[1])]
            continue
        if g=='ry':mat=np.array([[np.cos(a/2),-np.sin(a/2)],[np.sin(a/2),np.cos(a/2)]])
        elif g=='rz':mat=np.diag([np.exp(-1j*a/2),np.exp(1j*a/2)])
        else:mat=fixed[g]
        lo,hi=pairs[qs[0]];v,w=state[lo].copy(),state[hi].copy()
        state[lo]=mat[0,0]*v+mat[0,1]*w;state[hi]=mat[1,0]*v+mat[1,1]*w
    return state


def edges(L,m):
    # Independent layer generation, not exporter.bonds.
    if m[0]=='x':return [(y*L+x,y*L+x+1) for y in range(L) for x in range(int(m[-1]),L-1,2)]
    return [(y*L+x,(y+1)*L+x) for y in range(int(m[-1]),L-1,2) for x in range(L)]


def spin_state(L,rows):
    n=L*L;ids=np.arange(2**n);v=np.ones(2**n,complex)
    for i,j in edges(L,'x_0'):
        bi=(ids>>i)&1;bj=(ids>>j)&1
        v*=np.where(bi!=bj,1-2*bi,0)/np.sqrt(2)
    for row in rows:
        if row['gate']!='eSWAP':continue
        a=float(row['angle_rad'])/2
        for i,j in edges(L,row['matching']):
            diff=((ids>>i)^(ids>>j))&1
            perm=ids^(diff<<i)^(diff<<j)
            v=np.cos(a)*v-1j*np.sin(a)*v[perm]
    return v


def distance(a,b):
    phase=np.vdot(b,a)/np.vdot(b,b)
    assert abs(abs(phase)-1)<1e-10
    return float(np.max(np.abs(a-phase*b)))


def main():
    check_local()
    report={'checks':[], 'scope':'Logical, unsynthesized circuits; no large-Hubbard state-vector simulation.'}
    data={}
    for model in ('heisenberg','hubbard'):
        with (ROOT/f'{model}_layer_shared.csv').open() as f:data[model]=list(csv.DictReader(f))
        for nb in sorted({int(r['N_b']) for r in data[model]}):
            rows=[r for r in data[model] if int(r['N_b'])==nb]
            expected=[('eSWAP',b,m) for b in range(1,nb+1) for m in ORDER]
            if model=='hubbard':expected += [('SW',1,m) for m in ORDER]
            assert [(r['gate'],int(r['block']),r['matching']) for r in rows]==expected
    manifest=json.loads((ROOT/'qasm/manifest.json').read_text())
    for entry in manifest:
        path=ROOT/entry['file'];n,ops=parse(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']
        L,nb,model=entry['L'],entry['N_b'],entry['model'];hub=model=='hubbard'
        assert n==L*L*(2 if hub else 1)
        actual=[a for g,_,a in ops if g in ('ry','rz')]
        assert len(actual)==(3*nb+8*hub)*2*L*(L-1)==entry['rotations']
        rows=[r for r in data[model] if int(r['N_b'])==nb]
        expected=[]
        for r in rows:
            a=float(r['angle_rad'])
            angles=[a/2,-a/2,-a/2] if r['gate']=='eSWAP' else [a/8,-a/8,-a/8,a/8]*2
            expected.extend(angles*len(edges(L,r['matching'])))
        assert np.array_equal(actual,expected)
        if L==4:
            if hub:
                prefix=ops[:4*(L*L//2)+11*nb*2*L*(L-1)]
                assert all(q%2==0 for _,qs,_ in prefix for q in qs)
                got=evolve([(g,tuple(q//2 for q in qs),a) for g,qs,a in prefix],16)
                saved=np.load(ROOT/f'hubbard/L4_parameters/Nb{nb}/best.npz')['spin_statevector']
                err=distance(got,saved)
                assert err<1e-10,(model,nb,err)
            else:
                got=evolve(ops,16);err=distance(got,spin_state(4,rows))
                assert err<1e-10,(model,nb,err)
            report['checks'].append({'model':model,'L':L,'N_b':nb,'max_state_error_up_to_phase':err,
                                     'test':'saved pre-SW spin state' if hub else 'full exported circuit vs direct eSWAP evolution'})
    # Full Hubbard circuit including embedding and the complete SW block on 2x2.
    for nb in (3,4,5):
        rows=[r for r in data['hubbard'] if int(r['N_b'])==nb]
        got=evolve(circuit(2,rows,True),8)
        spin=spin_state(2,rows);target=np.zeros(256,complex)
        for bits,a in enumerate(spin):
            occupied=sum(1<<(2*i+((bits>>i)&1)) for i in range(4))
            target[occupied]=a
        for r in rows:
            if r['gate']=='SW':
                for i,j in edges(2,r['matching']):target=target_sw(i,j,8,float(r['angle_rad']))@target
        err=distance(got,target);assert err<1e-10,(nb,err)
        report['checks'].append({'model':'hubbard','L':2,'N_b':nb,'test':'full circuit vs fermionic SW evolution','max_state_error_up_to_phase':err})
    report['exported_files_checked']=len(manifest)
    (ROOT/'qasm/validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
