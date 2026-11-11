[English](Blast-Wave) · **日本語** · ソース: [`examples/blast_wave`](https://github.com/ytworks/pypic/tree/main/examples/blast_wave)

# 磁化天体のまわりの爆風

高密度の電子–イオンプラズマの殻（デブリ）が、開いた 3 次元の箱の中心から、静止した磁化プラズマ
（アルヴェン・マッハ数 {{ c.gate.params.mach_alfven }}、m_i/m_e = {{ c.gate.params.mass_ratio }}）の中へ膨張します。
殻は磁場と周囲のプラズマを外へ押し、後ろに磁気空洞を残します。中心には磁気双極子をもつ吸収球があり、殻は
その外側から始まります。この run では、圧縮された磁場の波面は速い磁気音波の速度の
{{ r.gate.front_mach_fast | .2f }} 倍で進み、デブリのイオンは磁場を横切る方向に
{{ r.gate.debris_extent | .1f }} d_e まで広がります。これは Le et al. (2021) のスケーリングの
{{ r.gate.debris_extent_ratio | .2f }} 倍です。

## 物理

**膨張のスケール**（Le et al. 2021）。質量 M_d のデブリの雲が、一様磁場 B₀ の中の密度 n₀ のプラズマへ、中心から
速度 V_d で広がります。等質量半径 R_m は、その内側の周囲プラズマの質量がデブリの質量と等しくなる半径です。
磁気的な停止半径 R_B は、その内側の磁場のエネルギーがデブリの運動エネルギーと等しくなる半径です。

```math
\frac{4\pi}{3} R_m^3\, n_0 m_i = M_d , \qquad
\frac{4\pi}{3} R_B^3\, \frac{B_0^2}{2} = \frac{1}{2} M_d V_d^2
\;\;\Rightarrow\;\; R_B = R_m M_A^{2/3} , \quad M_A = \frac{V_d}{v_A} .
```

デブリのイオンは、流れの速度でのジャイロ半径 ρ_d = V_d/Ω_ci が R_m より小さいとき、場を通して周囲のプラズマと
結合します。ρ_d > R_m では結合しません（Le et al. 2021）。

**デブリの広がり。** Le et al. (2021) は、磁場を横切る方向に線積分したデブリ密度の 98 % を含む半径 Δ98 を測り、
その最大値を次の式で近似しています。

```math
\frac{\Delta_{max}}{R_m} \sim \frac{1 + c_0 M_A}{\left(1 + c_1 / M_A^2\right)^{1/D}} ,
\qquad c_0 = 0.03 , \quad c_1 = 2 .
```

係数は 2 次元のハイブリッド計算（D = 2、ガウス型の雲、デブリと周囲のイオンの質量比 3）への当てはめです。
D = 3 について論文は式の形だけを示していて、当てはめはありません。下の比較は、この外挿との比較です。

**波面。** デブリの前方の圧縮は B₀ を横切って進みます。これを速い磁気音波の速度、および垂直 MHD 衝撃波の
圧縮比 X と比べます。

```math
v_f^2 = v_A^2 + \frac{\gamma (T_e + T_i)}{m_i} , \qquad
(2-\gamma) X^2 + \left[\gamma(\beta+1) + (\gamma-1) M^2\right] X - (\gamma+1) M^2 = 0 , \qquad \gamma = 5/3 .
```

M は波面のアルヴェン・マッハ数です。

| 量 | 値 |
|---|---|
| デブリの速度 V_d | {{ c.gate.planner.v_d | .4f }} c |
| アルヴェン速度、速い磁気音波の速度 | {{ c.gate.planner.v_a | .4f }} c、{{ c.gate.planner.v_fast | .4f }} c |
| プラズマ β（電子とイオン） | {{ c.gate.planner.beta | .2f }} |
| R_m、R_B | {{ c.gate.planner.r_m | .1f }} d_e、{{ c.gate.planner.r_b | .1f }} d_e |
| デブリのイオンのジャイロ半径 ρ_d | {{ c.gate.planner.debris_gyroradius | .1f }} d_e |
| D = 3 のスケーリングによる Δ_max | {{ c.gate.planner.debris_extent_le | .1f }} d_e |

## 実験条件

規格化単位: c = ε₀ = 1、長さは d_e = c/ω_pe、時間は 1/ω_pe、周囲の密度 n₀ = 1。B₀ と双極子モーメントは
+z を向きます。

- **4 つの種。** 周囲の電子とイオンは静止して箱を満たし、6 つの開放面から入ります。デブリの電子とイオンは、
  半径方向にガウス型の殻（中心 {{ c.gate.planner.shell_radius | .1f }} d_e、幅
  {{ c.gate.planner.shell_width | .1f }} d_e、最大密度 {{ c.gate.planner.shell_peak_density | .1f }} n₀）として
  置き、動径方向の速度 V_d と周囲と同じ温度を与えます。デブリのイオンの質量は周囲のイオンと同じです。
- **外部場**: B₀ と双極子は外部場です。格子はそれらとの差だけを持ち、出力は合計の場です。双極子の赤道面の
  磁場は B₀ と逆向きで、合計の磁場は中心から {{ c.gate.planner.dipole_radius | .1f }} d_e で 0 になります。
- **球**: 半径 {{ c.gate.planner.body_radius | .2f }} d_e。入った粒子を取り除き、その電荷を保持します
  （浮遊電位の絶縁体）。
- **run の長さ**: 1.25 R_B/V_d で、Ω_ci t = {{ c.gate.planner.omega_ci_t_end | .1f }} です。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![爆風](images/blast_wave/gate.gif)

球の内部は表示していません。

- **左、z = 中心の面の B_z / B₀**（B₀ に垂直な面）: 波面の圧縮された磁場と、その後ろの空洞。
- **中央、z = 中心の面の周囲イオンの密度。**
- **右、y = 中心の面のデブリのイオン密度**（B₀ を含む面、対数のカラースケール）。

## スケーリングとの比較

![波面、デブリ、空洞](images/blast_wave/gate-blast_wave.png)

B₀ に垂直な面の B_z と 8 本の線上の波面、B₀ を含む面のデブリのイオン、および波面・デブリ・空洞の半径の
時間変化です。波面は B_z > 1.25 B₀ となる最も外側の半径です。波面の基準には、波面が 0.3〜0.7 R_B にある
出力を使います。

| 量 | PIC | 参照値 | ゲート |
|---|---|---|---|
| 区間内の波面の速度 | {{ r.gate.front_speed | .4f }} c | {{ r.gate.fast_speed | .4f }} c（速い磁気音波の速度） | より速い |
| 区間内の B_z / B₀ の円周平均の最大値 | {{ r.gate.compression | .2f }} | — | > 1.5 |
| 波面の半径、軸方向と対角方向の差 | {{ r.gate.asymmetry | .2% }} | — | 5 % 以内 |
| デブリのイオンの Δ98 の最大値 | {{ r.gate.debris_extent | .1f }} d_e | {{ r.gate.debris_extent_le | .1f }} d_e（Le et al. 2021、D = 3） | 25 % 以内、終了より前 |
| 波面より外の n_e、n_i、\|B\| | ステップ 0 から {{ r.gate.ambient_deviation | .2% }} 以内 | — | 5 % |
| ガウスの法則、max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | — | ≤ 1e-10 |
| max\|∇·B\| | {{ r.gate.divB_max | .0e }} | — | ≤ 1e-12 |
| 測った波面の速度での垂直 MHD 衝撃波の圧縮比 | — | {{ r.gate.compression_mhd_front | .2f }} | 報告のみ |
| 磁気空洞（B_z < 0.5 B₀）の最大半径 | {{ r.gate.cavity_max | .2f }} d_e | — | 報告のみ |

ゲートの基準は最初の run の前に決めました。

## 補足

- **デブリの広がりはスケーリングより大きい値です。** スケーリングは 2 次元の当てはめを D = 3 で使ったもので、
  この例の雲は周囲と同じイオン質量の殻で、R_m は 8 d_i です。この run ではこれらの効果を分離できません。
- **波面はデブリから離れます。** run の終わりに波面は {{ r.gate.front_end | .1f }} d_e にあり、R_B より外です。
  デブリのイオンは、磁場を横切る方向には最大の広がりより内側にとどまります。B₀ に沿う方向にはより遠くまで
  広がり、中心からの rms 距離は B₀ に沿って {{ r.gate.debris_rms_along | .1f }} d_e、横切る方向に
  {{ r.gate.debris_rms_across | .1f }} d_e です。
- **圧縮された磁場は、波面が外へ進むにつれて弱くなります。** 表の値は区間内の最大値です。
- **エネルギー。** 終了時に、デブリのイオンは最初の運動エネルギーの
  {{ r.gate.debris_ion_energy_left | .0% }} を持ちます。周囲のイオンは、ステップ 0 のデブリ（電子とイオン）の
  運動エネルギーの {{ r.gate.ambient_ion_energy_gain | .0% }} を得ています。全体の収支は監視のみです
  （explicit スキーム）。変化は、同じデブリのエネルギーの {{ r.gate.energy_budget_change | +.1% }} です。
- **球**は終了までに {{ r.gate.n_absorbed }} 個の粒子を吸収し、電荷 {{ r.gate.body_charge | +.1f }}
  （e n₀ d_e³ 単位）を持ちます。デブリのない静止プラズマの試験では、電荷が一定値に落ち着くことを確認できません
  でした。球の電位は報告のみで、ゲートの基準ではありません。
- **粒子数**: 粒子配列は、最も埋まったランクで容量の {{ r.gate.fill_max | .0% }} まで埋まります。
- **run ごとの差**: この run と、同じシードの以前の FP64 の run 3 回で、波面の速度、円周平均、対称性は
  全桁で、Δ98 は 12 桁で一致しました。ほかの値の差は、値の 4e-5 以下でした。`precision: mixed` の run も
  同じ基準を満たしました。その記録は、リポジトリのゲート記録にあります。

## 再現方法

```bash
python -m pypic.launch -n 8 examples/blast_wave/run.py \
    --config examples/blast_wave/configs/gate.yaml --out output/wiki/blast_wave/gate
python examples/blast_wave/analyze.py output/wiki/blast_wave/gate --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- A. Le, D. Winske, A. Stanier, W. Daughton, M. Cowee, B. Wetherton and F. Guo, "Astrophysical
  explosions revisited: collisionless coupling of debris to magnetized plasma", J. Geophys. Res. Space
  Physics 126, e2021JA029125 (2021).
- B. D. Keenan, A. Le, D. Winske, A. Stanier, B. Wetherton, M. Cowee and F. Guo, "Hybrid
  particle-in-cell simulations of electromagnetic coupling and waves from streaming burst debris",
  Phys. Plasmas 29, 012107 (2022).
- M. E. Dieckmann et al., "Expansion of a radially symmetric blast shell into a uniformly magnetized
  plasma", Phys. Plasmas 25, 052108 (2018).
- H. B. Nersisyan and D. A. Osipyan, "Collisionless plasma expansion in the presence of a dipole
  magnetic field", Contrib. Plasma Phys. 49, 351 (2009).
