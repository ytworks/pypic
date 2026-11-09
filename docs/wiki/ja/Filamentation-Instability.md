[English](Filamentation-Instability) · **日本語** · ソース: [`examples/filamentation`](https://github.com/ytworks/pypic/tree/main/examples/filamentation)

# 電流フィラメンテーション（ワイベル）不安定性

±x 方向に互いを通り抜ける 2 本のプラズマビームは、流れに垂直な摂動に対して不安定です。電流の小さな横方向の
さざ波が磁場を作り、その磁場が同じ向きの電流を引き寄せ、逆向きの電流を押し離します。さざ波は電流
フィラメントに成長し、フィラメントの運ぶ磁場が、のちに無衝突衝撃波を媒介します。ペアプラズマで測った
モード m = 2…10 の成長率は、冷たい相対論的理論から中央値で {{ r.pair.median_rel_error_m2_10 | .1% }} の
ずれです（ゲート 10 %）。

## 物理

x 方向にドリフト β（ローレンツ因子 γ）で進む冷たいビームと、流れに垂直な波数ベクトル k に対して、相対論的
流体理論は純成長モードの成長率 Γ を与えます（c = 1）。

```math
\Gamma^4 + \Gamma^2\left(k^2 + \sum_s \frac{\omega_{ps}^2}{\gamma_s^3}\right) - k^2 \sum_s \frac{\omega_{ps}^2\beta_s^2}{\gamma_s} = 0 .
```

ドリフト方向の応答には縦方向の質量 γ³m が効きます。Γ は k とともに増え、k ≫ ω_p/c で βω_p/√γ に近づきます。
プランナーは計算領域に入る全フーリエモードについて Γ を計算します。

| 設定 | ビーム | 領域内の最大成長モードの Γ | k → ∞ の Γ | 捕捉による B_sat の見積もり |
|---|---|---|---|---|
| `pair` | 電子と陽電子がそれぞれ u = ±{{ c.pair.planner.u0 }} の 2 本のビーム | {{ c.pair.planner.gamma_max | .3f }} ω_pe | {{ c.pair.planner.gamma_limit | .3f }} ω_pe | {{ c.pair.planner.B_sat_trapping | .3f }} |
| `electron_ion` | 静止イオン上の u = ±{{ c.electron_ion.planner.u0 }} の電子ビーム | {{ c.electron_ion.planner.gamma_max | .3f }} ω_pe | {{ c.electron_ion.planner.gamma_limit | .3f }} ω_pe | {{ c.electron_ion.planner.B_sat_trapping | .3f }} |

磁場は、フィラメントの磁場がビーム粒子を捕捉するとき、つまり磁場中のバウンス周波数が成長率に達するときに
おおよそ飽和します: B_sat ≈ γ₀mΓ²/(qβ₀k)（Davidson et al. 1972）。

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe。計算領域は流れの方向に 4 セルしかなく、二流体モード
（k ∥ v）を除外しています。(y, z) 面内の k を、スキン長 c/ω_p あたり 10 セルで分解します。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![ペアプラズマのフィラメンテーション](images/filamentation/pair.gif)

`pair`、流れに垂直な (y, z) 面。

- **左上、B_z** と **右上、J_x**: 雑音から、符号の交互に並ぶ電流フィラメントが現れ、それぞれが自分の磁場に
  包まれます。その後、小さなフィラメントが合体して、数の少ない大きなフィラメントになります（飽和後の合体）。
- **左下、位相空間 (y, u_x)**: u_x = ±{{ c.pair.planner.u0 }} の 2 本のビームは、磁場の成長とともに曲げられます。
  後半では、粒子がフィラメントの磁場に捕捉されてフィラメントを横切るように振動します。
- **右下、エネルギー**: 磁場エネルギーは成長率の 2 倍で指数関数的に増え、初期運動エネルギーの約
  {{ r.pair.eps_B_peak | .0% }} で飽和します。

## 線形理論との比較

![k に対する成長率（ペアプラズマ）](images/filamentation/pair-filamentation.png)

`pair`。**左**: 共通の線形区間で |B(k_⊥)| から測った各モードの成長率（点）と、冷たい理論の曲線 Γ(k)（破線）。
**中央**: 3 つのモードの |B(k)| の時間変化。速いモードほど先に飽和します。**右**: 磁場エネルギーが最大の
ときの B_z。

![k に対する成長率（電子–イオン）](images/filamentation/electron_ion-filamentation.png)

`electron_ion`（報告のみで、ゲート設定ではありません）。同じ解析です。

| 設定 | \|Γ_PIC − Γ_理論\|/Γ_理論 の中央値（m = 2…10） | 最大 | ゲート | 線形区間 [1/ω_pe] | 最大時の ε_B | 最大の時刻 [1/ω_pe] |
|---|---|---|---|---|---|---|
| `pair` | {{ r.pair.median_rel_error_m2_10 | .1% }} | {{ r.pair.max_rel_error_m2_10 | .1% }} | 10 % | {{ r.pair.linear_window.0 | .1f }} – {{ r.pair.linear_window.1 | .1f }} | {{ r.pair.eps_B_peak | .3f }} | {{ r.pair.t_peak | .1f }} |
| `electron_ion` | {{ r.electron_ion.median_rel_error_m2_10 | .1% }} | {{ r.electron_ion.max_rel_error_m2_10 | .1% }} | 報告のみ | {{ r.electron_ion.linear_window.0 | .1f }} – {{ r.electron_ion.linear_window.1 | .1f }} | {{ r.electron_ion.eps_B_peak | .3f }} | {{ r.electron_ion.t_peak | .1f }} |

ε_B は磁場エネルギーを初期運動エネルギーで割ったものです。ガウスの法則は {{ r.pair.gauss_max | .1e }}
（`pair`）と {{ r.electron_ion.gauss_max | .1e }}（`electron_ion`）の精度で成り立っています。

## 補足

- **共通のフィット区間.** 多くのモードが同時に成長し、最も速いモードが飽和するとすべてのモードが線形領域を
  外れます。そこで区間は全磁場エネルギーから決めます。振幅で初期の雑音の e² 倍になった時刻から、ピークの
  1/e² になる時刻（ピークの手前）までです。
- **モード m = 1** はゲートから除外しています。最も成長が遅く、速いモードが飽和する時点でまだ雑音に
  埋もれているためです。
- **`electron_ion`** は報告のみです。測定した成長率は小さい k では理論曲線に沿いますが、大きい k では系統的に
  下回ります（ずれの中央値 {{ r.electron_ion.median_rel_error_m2_10 | .1% }}）。原因はまだ調べていません。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/filamentation/run.py --config examples/filamentation/configs/pair.yaml \
    --out output/wiki/filamentation/pair
python examples/filamentation/analyze.py output/wiki/filamentation/pair --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- E. S. Weibel, "Spontaneously growing transverse waves in a plasma due to an anisotropic velocity
  distribution", Phys. Rev. Lett. 2, 83 (1959).
- B. D. Fried, "Mechanism for instability of transverse plasma waves", Phys. Fluids 2, 337 (1959).
- R. C. Davidson, D. A. Hammer, I. Haber and C. E. Wagner, "Nonlinear development of electromagnetic
  instabilities in anisotropic plasmas", Phys. Fluids 15, 317 (1972).
