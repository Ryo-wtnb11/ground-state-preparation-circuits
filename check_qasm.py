"""Numerical checks for the export decomposition (requires NumPy)."""
import re
import numpy as np
def expm(a):
    # All test generators are anti-Hermitian.
    w,v=np.linalg.eigh(1j*a)
    return (v*np.exp(-1j*w))@v.conj().T

from export_qasm import ROOT, eswap, sw


def simulate(ops, n):
    u = np.eye(2**n, dtype=complex)
    for g,qs,a in ops:
        if g == 'cx':
            indices=np.arange(2**n)
            perm=indices ^ (((indices >> qs[0]) & 1) << qs[1])
            u=u[perm]
            continue
        matrices={'h':np.array([[1,1],[1,-1]])/np.sqrt(2),
                  'x':np.array([[0,1],[1,0]]), 'z':np.diag([1,-1]),
                  's':np.diag([1,1j]), 'sdg':np.diag([1,-1j])}
        if g=='rz': m=np.diag([np.exp(-1j*a/2),np.exp(1j*a/2)])
        elif g=='ry': m=np.array([[np.cos(a/2),-np.sin(a/2)],[np.sin(a/2),np.cos(a/2)]])
        else: m=matrices[g]
        for b in range(2**n):
            if b & (1<<qs[0]): continue
            c=b|(1<<qs[0]);u[[b,c]]=m@u[[b,c]]
    return u


def target_sw(i,j,n,phi):
    h=np.zeros((2**n,2**n)); doublons=[]
    for b in range(2**n):
        doublons.append(sum(((b>>(2*k))&1)*((b>>(2*k+1))&1) for k in (i,j)))
        for spin in (0,1):
            for dst,src in ((2*i+spin,2*j+spin),(2*j+spin,2*i+spin)):
                if not (b>>src)&1 or (b>>dst)&1: continue
                rem=b^(1<<src)
                parity=(b&((1<<src)-1)).bit_count()+(rem&((1<<dst)-1)).bit_count()
                h[rem|(1<<dst),b]+=(-1)**parity
    q=np.array(doublons)
    return expm(phi/4*(q[:,None]-q[None,:])*h)


def check():
    swap=np.eye(4)[[0,2,1,3]]
    for a in (0.,.427,-1.23):
        got=simulate(eswap(0,1,a),2); target=expm(-1j*a/2*swap)
        phase=np.vdot(target,got)/4
        assert np.allclose(got,phase*target,atol=1e-12)
        for j in (1,2):
            assert np.allclose(simulate(sw(0,j,a),2*(j+1)),target_sw(0,j,2*(j+1),a),atol=1e-12)
    count=0
    for path in sorted((ROOT/'qasm').rglob('*.qasm')):
        lines=path.read_text().splitlines(); assert lines[:2]==['OPENQASM 2.0;','include "qelib1.inc";']
        n=int(re.fullmatch(r'qreg q\[(\d+)\];',lines[4])[1])
        for line in lines[5:]:
            assert re.fullmatch(r'(h|x|z|s|sdg|cx|r[yz]\([-+0-9.e]+\)) q\[\d+\](,q\[\d+\])?;',line),line
            assert all(int(q)<n for q in re.findall(r'q\[(\d+)\]',line))
        count+=1
    print(f'PASS: eSWAP matrix up to global phase; SW including intermediate JW sites; {count} QASM files syntax and indices.')

if __name__=='__main__': check()
