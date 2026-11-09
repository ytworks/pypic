[English](Plasma-Expansion) · **日本語** · ソース: [`examples/plasma_expansion`](https://github.com/ytworks/pypic/tree/main/examples/plasma_expansion)

# 真空へのプラズマ膨張（Mora 2003）

熱い電子と冷たいイオンからなる板状のプラズマを、真空中へ自由に膨張させます。電子が先行してイオン前面に
電荷分離のシースを作り、その電場がイオンを加速します。これは高強度レーザーによるイオン加速の機構である
TNSA（target-normal sheath acceleration）です。Mora（2003）は等温電子の場合の膨張を厳密に解きました。
陽子質量の場合、PIC のイオンは予言どおりの自己相似な速度分布に
{{ r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }} の精度で従い、最大エネルギーは等温則より
{{ r.isothermal_hydrogen.E_max_median_rel_error | .0% }} 低く出ます。PIC の運動論的な電子は、モデルの理想的な
電子より冷たいためです。

## 物理

τ = ω_pi t/√(2e)、イオン価数 Z、電子温度 T_e、初期デバイ長 λ_D0 とすると、等温の自己相似解は次を与えます。

```math
E_{max} = 2 Z T_e \left[\ln\left(\tau + \sqrt{\tau^2 + 1}\right)\right]^2 ,
\qquad
x_f = 2\sqrt{2e}\,\lambda_{D0}\left[\tau \ln\left(\tau + \sqrt{\tau^2+1}\right) - \sqrt{\tau^2+1} + 1\right],
```

```math
\frac{dN}{dE} = \frac{n_{i0}\, c_s t}{\sqrt{2 E Z T_e}} \exp\left(-\sqrt{\frac{2E}{Z T_e}}\right),
\qquad v = c_s + \frac{x}{t}\ \ \text{(quasi-neutral rarefaction)} ,
```

ここで c_s = √(ZT_e/m_i) はイオン音速です。前面の電場は E_f = 2E₀/√(2e + ω_pi²t²) で減衰します。

| 設定 | m_i/m_e | 板の厚さ | 計算終了 | 用途 |
|---|---|---|---|---|
| `isothermal_hydrogen` | {{ c.isothermal_hydrogen.params.mass_ratio | .0f }} | {{ c.isothermal_hydrogen.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal_hydrogen.params.tau_end }} | ゲート設定 |
| `isothermal` | {{ c.isothermal.params.mass_ratio }} | {{ c.isothermal.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal.params.tau_end }} | 高速。前面の速度が電子の熱速度に近づく |
| `thin_foil` | {{ c.thin_foil.params.mass_ratio }} | {{ c.thin_foil.params.slab_lambda_D }} λ_D | τ = {{ c.thin_foil.params.tau_end }} | 薄膜の膨張で電子が冷える。等温則は上限になる |
| `isothermal_hydrogen_long` | {{ c.isothermal_hydrogen_long.params.mass_ratio | .0f }} | {{ c.isothermal_hydrogen_long.params.slab_lambda_D }} λ_D | τ = {{ c.isothermal_hydrogen_long.params.tau_end }} | 記録済みの本番実行 |

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、エネルギーは m_e c²。T_e =
{{ c.isothermal_hydrogen.params.Te }} m_e c²（5.11 keV）なので λ_D0 = {{ c.isothermal_hydrogen.planner.lambda_D0 | .2f }} c/ω_pe です。
格子は準 1 次元（横方向 4 × 4 セル）で Δx = λ_D0、電子とイオンは同じ階段状の密度分布を共有するので、初期に
プラズマは厳密に中性です。x 方向の両端の壁は粒子を反射します。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![プラズマ膨張](images/plasma_expansion/isothermal_hydrogen_long.gif)

`isothermal_hydrogen_long`、τ = 30 までの記録済み本番実行。

- **左上、イオンの位相空間 (x, u_x)**: 板の端から希薄波が内側へ進み、端より外のイオンは直線状のランプ、
  つまり自己相似な速度分布 v = c_s + x/t に加速されます。
- **右上、イオンと電子の密度**: 密度は真空に向かって滑らかに減ります。電子とイオンは最前面を除いてほぼ重なって
  います（準中性）。
- **左下、E_x**: シースの電場はイオン前面でピークをもち、前面が外へ進むにつれて減衰します。
- **右下、イオンの最大エネルギー**: PIC で測った E_max(τ)（点）を、時間とともに Mora の等温則（破線）の上に
  描きます。

## Mora の解との比較

![膨張（陽子質量）](images/plasma_expansion/isothermal_hydrogen-expansion.png)

`isothermal_hydrogen`（ゲート）。**左**: E_max(τ)。**中央**: イオン前面の位置。**右**: 終了時のイオンの
エネルギースペクトルと Mora の dN/dE。

![膨張（長時間の run）](images/plasma_expansion/isothermal_hydrogen_long-expansion.png)

`isothermal_hydrogen_long`、τ = 30 まで。スペクトルはエネルギーで 3 桁にわたって自己相似な dN/dE に従います。

| 設定 | τ_end での E_max（PIC / Mora）[m_e c²] | E_max 誤差の中央値 | 前面位置の誤差の中央値 | 速度分布の誤差の中央値 | T_e 終了時 / 初期 | ゲート |
|---|---|---|---|---|---|---|
| `isothermal_hydrogen` | {{ r.isothermal_hydrogen.E_max_final | .4f }} / {{ r.isothermal_hydrogen.E_max_theory_final | .4f }} | {{ r.isothermal_hydrogen.E_max_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.front_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal_hydrogen.electron_temperature_final_over_initial | .2f }} | 速度 ≤ 10 %、E_max ≤ 20 % |
| `isothermal` | {{ r.isothermal.E_max_final | .4f }} / {{ r.isothermal.E_max_theory_final | .4f }} | {{ r.isothermal.E_max_median_rel_error | .1% }} | {{ r.isothermal.front_median_rel_error | .1% }} | {{ r.isothermal.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal.electron_temperature_final_over_initial | .2f }} | 報告のみ |
| `isothermal_hydrogen_long` | {{ r.isothermal_hydrogen_long.E_max_final | .4f }} / {{ r.isothermal_hydrogen_long.E_max_theory_final | .4f }} | {{ r.isothermal_hydrogen_long.E_max_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.front_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.velocity_profile_median_rel_error | .1% }} | {{ r.isothermal_hydrogen_long.electron_temperature_final_over_initial | .2f }} | 報告のみ |
| `thin_foil` | {{ r.thin_foil.E_max_final | .4f }} / {{ r.thin_foil.E_max_theory_final | .4f }} | {{ r.thin_foil.E_max_median_rel_error | .1% }} | {{ r.thin_foil.front_median_rel_error | .1% }} | {{ r.thin_foil.velocity_profile_median_rel_error | .1% }} | {{ r.thin_foil.electron_temperature_final_over_initial | .2f }} | E_max ≤ 等温則の 1.05 倍 |

誤差は 2 ≤ τ ≤ τ_end での |PIC − Mora|/Mora の中央値です。厚い板の run ではどれも PIC の値が等温モデルより
低くなります。ガウスの法則は {{ r.isothermal_hydrogen.gauss_max | .1e }}（`isothermal_hydrogen`）の精度で
成り立っています。

![薄膜](images/plasma_expansion/thin_foil-expansion.png)

`thin_foil`: 板の厚さが {{ c.thin_foil.params.slab_lambda_D }} λ_D しかないので、電子はエネルギーの大部分をイオンに
渡し（T_e は初期値の {{ r.thin_foil.electron_temperature_final_over_initial | .0% }} まで下がる）、E_max は等温則を
大きく下回ります。断熱膨張で期待されるとおりです（Mora 2005）。

<details><summary>m_i/m_e = 100</summary>

![膨張（m_i/m_e = 100）](images/plasma_expansion/isothermal-expansion.png)

</details>

## PIC が等温モデルより低くなる理由

- モデルは、無限の熱浴につながった等温（ボルツマン）電子を仮定しています。無衝突の PIC 電子にはそのような
  熱浴がなく、膨張領域は板の中より冷たくなります。
- すべてのずれ（E_max、前面位置、速度分布）は、膨張領域の実効的な T_e が約 15 % 低いことと整合します。
- このずれは Δx = λ_D/2、粒子数 2 倍、真空領域 2.5 倍でも変わらないので（感度試験は
  [example の README](https://github.com/ytworks/pypic/tree/main/examples/plasma_expansion) に記載）、数値的な
  ものではありません。τ = 10 から τ = 30 まで同じ程度にとどまり、m_i/m_e = 100 では大きくなります。そこでは
  前面の速度が電子の熱速度に近づき、モデルの前提 v_ion ≪ v_te が成り立たなくなるためです。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/plasma_expansion/run.py \
    --config examples/plasma_expansion/configs/isothermal_hydrogen.yaml --out output/wiki/plasma_expansion/isothermal_hydrogen
python examples/plasma_expansion/analyze.py output/wiki/plasma_expansion/isothermal_hydrogen --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- P. Mora, "Plasma expansion into a vacuum", Phys. Rev. Lett. 90, 185002 (2003).
- P. Mora, "Thin-foil expansion into a vacuum", Phys. Rev. E 72, 056401 (2005).
