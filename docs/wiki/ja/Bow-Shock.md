[English](Bow-Shock) · **日本語** · ソース: [`examples/bow_shock`](https://github.com/ytworks/pypic/tree/main/examples/bow_shock)

# 磁化惑星のバウショック

磁化した電子–イオンの風（アルヴェン・マッハ数 {{ c.gate.params.mach_alfven }}、m_i/m_e =
{{ c.gate.params.mass_ratio }}）が開いた 3 次元の箱を流れ、中心に磁気双極子をもつ吸収球に当たります。
双極子の磁場は磁気圏界面で風を止めます。風は速い磁気音波より速いので、磁気圏界面の上流にバウショックが
立ち、その間に圧縮されたプラズマの層（シース）ができます。この run では、磁気圏界面は圧力平衡の距離の
{{ r.gate.magnetopause | .2f }} 倍、衝撃波は {{ r.gate.shock | .2f }} 倍の位置にあり、Spreiter et al. (1966) の
気体力学的な離れ距離は {{ r.gate.shock_spreiter | .2f }} 倍です。

## 物理

**磁気圏界面。** 風は、動圧が圧縮された双極子磁場の磁気圧と等しくなる位置で止まります（Chapman–Ferraro;
Spreiter et al. 1966）。平面の電流層の前面の磁場は双極子磁場 B_d = M/R³ の 2 倍なので、

```math
\frac{(2 B_d)^2}{2} = n\, m_i v_{sw}^2 , \qquad B_d = \frac{M}{R_{mp}^3} .
```

プランナーは距離 R_mp を受け取り、モーメント M を決めます。

**バウショック。** 惑星から風上へ向かう線（日下点線）の上では、惑星間磁場（IMF）は流れに垂直です。垂直 MHD
衝撃波のランキン–ユゴニオ条件から、密度の圧縮比 X = n₂/n₁ は次の式の正の根です。

```math
(2-\gamma) X^2 + \left[\gamma(\beta+1) + (\gamma-1) M_A^2\right] X - (\gamma+1) M_A^2 = 0 , \qquad \gamma = 5/3 ,
```

気体力学のモデルは、シースの厚さをこの圧縮比と結びつけます（Spreiter et al. 1966）。

```math
\frac{R_{bs} - R_{mp}}{R_{mp}} = \frac{1.1}{X} .
```

| 量 | 値 |
|---|---|
| プラズマ β（電子とイオン） | {{ c.gate.planner.beta | .2f }} |
| 速い磁気音波のマッハ数 | {{ c.gate.planner.mach_fast | .2f }} |
| MHD 衝撃波の圧縮比 X | {{ c.gate.planner.compression_mhd | .2f }} |
| シースの厚さ 1.1 / X | {{ c.gate.planner.sheath_spreiter | .2f }} R_mp |

**解析的でない既知の結果**: 衝撃波ができるには、磁気圏がイオンのスケールより大きいこと
（Dp = R_mp/d_i ≳ 1）が必要で、Dp ≈ 20 から地球型になります（Omidi et al. 2004）。Lavorenti et al. (2022) の
水星の PIC モデルは Dp ≈ 8 で、この run は Dp = {{ c.gate.params.dp }} です。

## 実験条件

規格化単位: c = ε₀ = 1、長さは d_e = c/ω_pe、時間は 1/ω_pe、風の密度 n₀ = 1。風は +x 方向に
v_sw = {{ c.gate.planner.v_sw | .3f }} c で進みます。IMF は +z、双極子モーメントは −z を向くので、昼側では
二つの磁場が平行です。R_mp = {{ c.gate.planner.r_mp | .0f }} d_e、球の半径は
{{ c.gate.planner.body_radius | .1f }} d_e です。

- **外部場**: IMF、運動電場 E = −v × B、双極子は外部場です。格子はそれらとの差だけを持つので、開放面の
  吸収層は場を風の値へ緩和します。出力は合計の場です。
- **開放面**: 6 面すべてが開放です。風は上流面（電子とイオンを同じ位置に対で入れる）と 4 つの側面から入り、
  下流面は粒子を取り除くだけです。
- **球**: 入った粒子を取り除き、その電荷を保持します（浮遊電位の絶縁体）。箱は最初、球の外が風で満たされて
  いて、双極子は t = 0 からあります。
- **run の長さ**: 通過時間 L_x/v_sw の {{ c.gate.params.crossings }} 倍で、IMF のイオン旋回周波数で
  Ω_ci t = {{ c.gate.planner.omega_ci_t_end | .1f }} です。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![バウショック](images/bow_shock/gate.gif)

風は左から入ります。球の内部は表示していません。

- **左、z = 中心の面のイオン密度**（IMF と双極子に垂直な面）。
- **中央、y = 中心の面のイオン密度**（IMF と双極子の軸を含む面）。
- **右、z = 中心の面の B_z / B_imf**（対数のカラースケール）。

## 圧力平衡・MHD 衝撃波との比較

![日下点線](images/bow_shock/gate-bow_shock.png)

二つの面のイオン密度（最後の半通過時間の平均）と、日下点線上の密度と磁気圧の分布です。磁気圏界面は、
密度が最大値の後ろで風の密度まで下がる位置、衝撃波は、密度が最初に風の密度の 2 倍に達する位置です。

![離れ距離](images/bow_shock/standoff.png)

| 量 | PIC | 参照値 | ゲート |
|---|---|---|---|
| 磁気圏界面の距離 / R_mp | {{ r.gate.magnetopause | .3f }} | 1（圧力平衡） | 25 % 以内 |
| 衝撃波の距離 / R_mp | {{ r.gate.shock | .3f }} | {{ r.gate.shock_spreiter | .3f }}（Spreiter et al. 1966） | 磁気圏界面より上流 |
| 日下点線上の最大密度 / n₀ | {{ r.gate.compression | .2f }} | {{ r.gate.compression_mhd | .2f }}（垂直 MHD 衝撃波） | ≥ 2 |
| シースの厚さ / 磁気圏界面の距離 | {{ r.gate.sheath | .2f }} | {{ r.gate.sheath_spreiter | .2f }}（Spreiter et al. 1966） | 報告のみ |
| 磁気圧で決めた磁気圏界面 / R_mp | {{ r.gate.magnetopause_magnetic | .3f }} | 1 | 報告のみ |
| 上流の密度と u_x（最後の 1 通過時間） | 風の値から {{ r.gate.upstream_deviation | .2% }} 以内 | — | 5 % |
| ガウスの法則、max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | — | ≤ 1e-10 |

ゲートの基準は最初の run の前に決めました。

## 補足

- **密度の最大値はシースの平均密度ではありません。** MHD の圧縮比より大きい値です。超臨界の垂直衝撃波は
  ランプの後ろでオーバーシュートすることが知られていて（Leroy et al. 1982）、この run のシースは厚さが
  約 2 d_i で、オーバーシュートの後ろに平坦部ができるには薄すぎます。この run では二つを分離できません。
- **シースは気体力学の値より厚い**: 磁気圏界面が R_mp の内側にあり、衝撃波は予測に近い距離にあるためです。
- **球は正に帯電します**: 吸収した電荷は終了時に {{ r.gate.body_charge | +.0f }}（e n₀ d_e³ 単位）です。吸収率は
  最初の通過時間の 1 ステップあたり {{ r.gate.absorbed_per_step.0 | .0f }} 個から 3 回目の
  {{ r.gate.absorbed_per_step.2 | .0f }} 個まで下がり、run の終わりでも一定になっていません。
- **粒子数**: 粒子配列（開始時の粒子数の 2 倍の容量）は、最も埋まったランクで容量の
  {{ r.gate.fill_max | .0% }} まで埋まります。
- **run ごとの差**: GPU の電流堆積はアトミック加算で、加算の順序が run ごとに変わります。同じシードで
  この設定を繰り返し実行すると、上の距離は 3 桁目が、球の電荷は 2 % 異なります（4 回の実行で +730 から +749）。
- **より大きい磁気圏**（`dp10`、記録つきの run）: R_mp = {{ c.dp10.params.dp }} d_i（上の run は
  {{ c.gate.params.dp }} d_i）、開始時の粒子数は {{ c.dp10.particles | , }} 個です。磁気圏界面は
  {{ r.dp10.magnetopause | .3f }} R_mp、衝撃波は {{ r.dp10.shock | .3f }} R_mp にあります（Spreiter et al.
  1966 の値は {{ r.dp10.shock_spreiter | .3f }}）。日下点線上の最大密度は {{ r.dp10.compression | .2f }} n₀、
  シースの厚さは磁気圏界面の距離の {{ r.dp10.sheath | .2f }} 倍です。この run はゲートの基準を満たします。
  stand-off 距離の図の 2 つ目の点がこの run です。

  <details><summary>`dp10` の日下点線</summary>

  ![dp10](images/bow_shock/dp10-bow_shock.png)

  </details>
- **エネルギー**は監視のみです（explicit スキーム）。

## 再現方法

```bash
python -m pypic.launch -n 8 examples/bow_shock/run.py \
    --config examples/bow_shock/configs/gate.yaml --out output/wiki/bow_shock/gate
python examples/bow_shock/analyze.py output/wiki/bow_shock/gate --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- J. R. Spreiter, A. L. Summers and A. Y. Alksne, "Hydromagnetic flow around the magnetosphere",
  Planet. Space Sci. 14, 223 (1966).
- N. Omidi, X. Blanco-Cano, C. T. Russell and H. Karimabadi, "Dipolar magnetospheres and their
  characterization as a function of magnetic moment", Adv. Space Res. 33, 1996 (2004).
- F. Lavorenti et al., "Electron dynamics in small magnetospheres: insights from global, fully kinetic
  plasma simulations of the planet Mercury", A&A 664, A133 (2022).
- M. M. Leroy, D. Winske, C. C. Goodrich, C. S. Wu and K. Papadopoulos, "The structure of perpendicular
  bow shocks", J. Geophys. Res. 87, 5081 (1982).
- E. R. Priest, *Solar Magnetohydrodynamics* (Reidel, 1982).
