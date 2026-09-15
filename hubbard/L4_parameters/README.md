# Hubbard L=4：最適化済みパラメーター

4×4 OBC（16サイト、32量子ビット）、t=1、U=8、Nup=Ndown=8。
N_b=3,4,5 の結果を2026-09-15に再計算（lg ジョブ 264568.CPUcluster-PRISM）。再対角化はしていない。

**Heisenberg回路のthetaは `heisenberg/parameters/l4_P{3,4,5}_layer_shared.npz`（= `heisenberg_layer_shared.csv`、論文 Table II と同一）に固定し、追加するSW回路の4つのphiのみを、Hubbard厳密基底状態へのFidelity最大化で最適化した結果。**
Hubbardエネルギーの最小化ではない。Fidelityは規格化した内積の絶対値二乗。
全ケースでBFGSの収束条件を満たしているが、大域的最適性を保証するものではない。

2026-09-08 の旧結果（ジョブ 264550）は N_b=4,5 の Heisenberg theta が別の最適化 run（Heisenberg Fidelity 0.9635 / 0.9859）に基づいていたため差し替えた。

## ファイル

- `summary.csv`：各N_bについてSWなし、全phi=0.5、最適化後のFidelityとエネルギー。
- `sw_phi.csv`：各N_bの最適SW角度（ラジアン）、列順yo,ye,xo,xe。
- `heisenberg_theta.csv`：固定Heisenberg角度（ラジアン）。blockは1始まり、各行の層順yo,ye,xo,xe。
- `parameters.json`：模型・回路条件、全thetaとphi、元データの場所を含む機械可読データ。
- `Nb3/`, `Nb4/`, `Nb5/`：元のbest.npz、config.json、result.json、evaluations.jsonl、energy.json。
- `sha256.json`：コピーした元ファイルのSHA256。

各`best.npz`はphi, theta, spin_statevector, fidelity, gradient, metadata_jsonを含む。
spin_statevectorはSW適用前の固定Heisenberg状態（65536成分）。
SW後の巨大なHubbard状態自体は保存していないが、この状態とphiで再構成できる。

## 結果

| N_b | SWなしFidelity | SW全角度0.5のFidelity | 最適化後Fidelity | 最適化後Hubbardエネルギー |
| --- | ---: | ---: | ---: | ---: |
| 3 | 0.42115357 | 0.83430031 | 0.85147944 | -6.34471356 |
| 4 | 0.44190426 | 0.87628921 | 0.89406118 | -6.36324125 |
| 5 | 0.45220662 | 0.89713327 | 0.91521018 | -6.37707549 |

エネルギーはH=-Σt(c†_iσ c_jσ+h.c.)+UΣn_i↑n_i↓の全系の値、定数シフトなし。
基底エネルギーE0=-6.808414466441853。
丸め前の値はsummary.csvを参照。

## 回路規約

- サイトi=y*Lx+x。軌道はup0,down0,up1,down1,...、qubit 2i=up、2i+1=down。
- Heisenberg spin bit 0=up、1=down。初期状態はxeのsinglet covering。
- 各Heisenbergブロックはyo→ye→xo→xe。thetaはブロックごとに4成分、各層内の結合で共有。
- theta配列の添字は0始まりで `theta[4*block+layer]`。
- SWはyo→ye→xo→xeの1 sweep。各層内でphiを共有。
- W_ij(phi)=exp(phi[Q_ij,h_ij]/4)=exp(-i phi omega_ij/2)。
  Q_ijは結合上のdoublon数、h_ijは正符号のhopping演算子、omega_ij=i[Q_ij,h_ij]/2。
- phi=4t/U=0.5が未最適化値。層間は可換とはしていない。
- SWの固定粒子数セクター実装はQulacsとの状態・解析勾配の一致を小系で検証済み。

## 出所

元の計算プロジェクト：
`/Users/ryowatanabe/Research/codes/StatePreparation/diagonalization/OBC/`

実装：hubbard_fidelity.py、optimize_sw_fidelity.py、sector_sw.py、sector_sw.cpp、prepare_fixed_heisenberg.py、convert_layer_shared.py。
ジョブ：264568（2026-09-15）、スクリプト `jobs/run_sw_layer_shared_large.sh`。
入力角度：`data/obc_Lx=4_Ly=4_P={3,4,5}_layer_shared.npz`（convert_layer_shared.py で本リポジトリの npz から変換）。
保存済みEDはlg上の
`/home/ryo-w/codes/StatePreparation/diagonalization/OBC/data/hubbard_4x4_U8_t1_tp1_sz0_Aplusplus_xdiag/`。
EDの巨大ファイルは本フォルダには複製していない。
