[English](Relativistic-Shock) · **日本語** · ソース: [`examples/relativistic_shock`](https://github.com/ytworks/pypic/tree/main/examples/relativistic_shock)

# 相対論的電子陽電子衝撃波

ローレンツ因子 γ₀ = {{ c.pair_2d.params.gamma0 }} の冷たい電子陽電子の流れを反射壁にぶつけます（ピストン法）。
入射するプラズマと反射したプラズマが互いを通り抜け、フィラメンテーション（ワイベル）不安定性が磁場を作り、
無衝突衝撃波ができて流れの上流へ進みます。衝撃波の内部で何が起きても、衝撃波をはさんだ速度と密度の跳びは
粒子数・運動量・エネルギーの保存から決まります。PIC の衝撃波は、2 次元でも 3 次元でも相対論的なジャンプ条件を
±10 % のゲート内で再現し、長時間の run ではフェルミ加速による非熱的なテイルができ始めます。

## 物理

壁の静止系が下流の静止系です。下流のプラズマを温度 θ₂ = T₂/(mc²) の相対論的マクスウェル–ユトナー（Synge）
気体とすると、エネルギーと粒子数の保存から次が得られます。

```math
\langle\gamma\rangle(\theta_2) = \gamma_0, \qquad
\beta_{sh} = \frac{\theta_2}{\gamma_0\beta_0}, \qquad
\frac{n_2}{n_1} = 1 + \frac{\beta_0}{\beta_{sh}} ,
```

⟨γ⟩(θ) は気体の平均ローレンツ因子で、粒子運動の自由度の数によって変わります。準 2 次元の run では 2、3 次元の
run では 3 です。

| 気体 | θ₂ | β_sh | n₂/n₁ |
|---|---|---|---|
| 2 次元 Synge 気体 | {{ c.pair_2d.planner.jump_2d.theta2 | .2f }} | {{ c.pair_2d.planner.jump_2d.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_2d.compression | .3f }} |
| 3 次元 Synge 気体 | {{ c.pair_2d.planner.jump_3d.theta2 | .2f }} | {{ c.pair_2d.planner.jump_3d.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_3d.compression | .3f }} |
| 断熱指数一定 Γ = 3/2 | — | {{ c.pair_2d.planner.jump_gamma_3_2.beta_sh | .3f }} | {{ c.pair_2d.planner.jump_gamma_3_2.compression | .3f }} |

**解析的でない既知の結果**: 衝撃波を何度も横切る粒子はエネルギーを得て（フェルミ加速）、べき乗則のテイル
dN/dγ ∝ γ^−p を作ります。2 次元では p ≈ 2.4（Spitkovsky 2008; Sironi et al. 2013）、等方拡散の理論では
p ≈ 2.2〜2.3（Kirk et al. 2000; Achterberg et al. 2001）で、最大エネルギーは γ_max ∝ t^{1/2} で増えます。

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe。相対論的スキン長は √γ₀ c/ω_pe =
{{ c.pair_2d.planner.skin | .3f }} です。流れは右から −x 方向に入り、壁は x = 0 にあります。補償付き二項電流
フィルタで、相対論的な流れの数値チェレンコフ不安定性を抑えます。有効な時間は、反対側の壁に開く真空の隙間で
t < L_x/(1 + β_sh) に制限されます。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![相対論的衝撃波](images/relativistic_shock/pair_2d.gif)

`pair_2d`、x–y 面（壁は左側）。

- **左上、B_z**: 反射した粒子が流れと出会う場所で磁場のフィラメントができ、衝撃波の後ろでは磁場が乱流状に
  なります。
- **右上、電子密度**: 圧縮された下流領域が壁から広がり、その端（衝撃波）が β_sh ≈ 0.5 で右へ進みます。
- **左下、電子の位相空間 (x, u_x)**: u_x = −{{ c.pair_2d.planner.u0 | .1f }} の入射ビームは衝撃波で熱化され、反射した
  粒子は正の u_x で衝撃波の前方を走ります。
- **右下、n_e(x, t)**: 原点から伸びる直線が衝撃波で、その傾きが衝撃波の速度です。

## ジャンプ条件との比較

### 全 run の衝撃波速度と圧縮比

![全 run のジャンプ条件](images/relativistic_shock/jump_conditions.png)

両方の量について、各 run の次元に合わせた気体（準 2 次元の run は 2 次元、`pair_3d` は 3 次元）の理論値に対する
測定値の比を示します。灰色の帯が ±10 % のゲートです。

| run | 気体 | β_sh PIC / 理論 | 誤差 | n₂/n₁ PIC / 理論 | 誤差 | ゲート |
|---|---|---|---|---|---|---|
| `pair_2d` | {{ r.pair_2d.dims_eff }} 次元 | {{ r.pair_2d.beta_sh | .3f }} / {{ r.pair_2d.beta_sh_theory | .3f }} | {{ r.pair_2d.beta_sh_rel_error | +.1% }} | {{ r.pair_2d.compression | .2f }} / {{ r.pair_2d.compression_theory | .2f }} | {{ r.pair_2d.compression_rel_error | +.1% }} | ±10 % |
| `pair_2d_long` | {{ r.pair_2d_long.dims_eff }} 次元 | {{ r.pair_2d_long.beta_sh | .3f }} / {{ r.pair_2d_long.beta_sh_theory | .3f }} | {{ r.pair_2d_long.beta_sh_rel_error | +.1% }} | {{ r.pair_2d_long.compression | .2f }} / {{ r.pair_2d_long.compression_theory | .2f }} | {{ r.pair_2d_long.compression_rel_error | +.1% }} | 報告のみ |
| `pair_2d_long_run2` | {{ r.pair_2d_long_run2.dims_eff }} 次元 | {{ r.pair_2d_long_run2.beta_sh | .3f }} / {{ r.pair_2d_long_run2.beta_sh_theory | .3f }} | {{ r.pair_2d_long_run2.beta_sh_rel_error | +.1% }} | {{ r.pair_2d_long_run2.compression | .2f }} / {{ r.pair_2d_long_run2.compression_theory | .2f }} | {{ r.pair_2d_long_run2.compression_rel_error | +.1% }} | 報告のみ |
| `pair_3d` | {{ r.pair_3d.dims_eff }} 次元 | {{ r.pair_3d.beta_sh | .3f }} / {{ r.pair_3d.beta_sh_theory | .3f }} | {{ r.pair_3d.beta_sh_rel_error | +.1% }} | {{ r.pair_3d.compression | .2f }} / {{ r.pair_3d.compression_theory | .2f }} | {{ r.pair_3d.compression_rel_error | +.1% }} | 報告のみ |

3 次元の run は 3 次元気体の値（β_sh ≈ 1/3、n₂/n₁ ≈ 4）に、準 2 次元の run は 2 次元気体の値（β_sh ≈ 1/2、
n₂/n₁ ≈ 3）に跳びます。衝撃波は、粒子がいくつの方向に動けるかを知っているのです。

### 衝撃波の構造と下流のスペクトル

![pair_2d_long](images/relativistic_shock/pair_2d_long-shock.png)

`pair_2d_long`（t = 2000）。**左**: n_e(x, t) と、ジャンプ条件から求まる衝撃波の軌跡（破線）。**中央**: 測定
時刻の密度分布と n₂/n₁ の水準。**右**: 下流のスペクトル。非熱的なテイルをもつ相対論的マクスウェル分布と、
べき乗則のフィット。

![pair_3d](images/relativistic_shock/pair_3d-shock.png)

`pair_3d`（3 次元、t = 600）。

<details><summary>pair_2d（ゲート）と長時間 run の 2 つ目の実現</summary>

![pair_2d](images/relativistic_shock/pair_2d-shock.png)
![pair_2d_long_run2](images/relativistic_shock/pair_2d_long_run2-shock.png)

</details>

### 粒子加速（報告のみ）

| run | 非熱的な割合 | テイルの指数 p（区間フィット） | γ_max ∝ t^α | 終了時の ε_B |
|---|---|---|---|---|
| `pair_2d_long` | {{ r.pair_2d_long.nonthermal_fraction | .1% }} | {{ r.pair_2d_long.tail_index_window_p | .2f }} ± {{ r.pair_2d_long.tail_index_window_err | .2f }} | α = {{ r.pair_2d_long.gamma_max_growth_exponent | .2f }} | {{ r.pair_2d_long.eps_B_end | .4f }} |
| `pair_2d_long_run2` | {{ r.pair_2d_long_run2.nonthermal_fraction | .1% }} | {{ r.pair_2d_long_run2.tail_index_window_p | .2f }} ± {{ r.pair_2d_long_run2.tail_index_window_err | .2f }} | α = {{ r.pair_2d_long_run2.gamma_max_growth_exponent | .2f }} | {{ r.pair_2d_long_run2.eps_B_end | .4f }} |
| 文献（2 次元、収束後） | — | ≈ 2.4 | α = 0.5 | — |

テイルはできていますが、t = 2000 ω_pe⁻¹（相対論的プラズマ時間 √γ₀/ω_pe の約 520 倍）ではまだ若く、その指数
（{{ r.pair_2d_long.tail_index_window_p | .1f }}）は収束値 2.4 よりずっと急です。Spitkovsky はこれを約 8700
ω_pe⁻¹ で測っています。収束したテイルには t ≳ 8000 ω_pe⁻¹ と、準 2 次元で約 10⁹ 個の粒子が必要です
（ロードマップ項目 R11）。

## 補足

- **長時間 run の 2 つの実現**: 同じ設定・同じシードで、この設定に影響しない 2 つのコミットで実行したもの
  です。GPU の電流割り付けは原子的加算の順序が非決定的で、乱流が丸め誤差の違いを増幅するため、2 つの run は
  独立な実現になります。その差は統計的な不確かさの目安です。
- **衝撃波速度のフィット区間**: 有効時間の 0.3〜0.7 倍。それより後では、反射粒子の前駆体が後退するプラズマの
  端に届いてしまいます。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/relativistic_shock/run.py \
    --config examples/relativistic_shock/configs/pair_2d.yaml --out output/wiki/relativistic_shock/pair_2d
python examples/relativistic_shock/analyze.py output/wiki/relativistic_shock/pair_2d --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- A. Spitkovsky, "On the structure of relativistic collisionless shocks in electron-ion plasmas",
  ApJ 673, L39 (2008); "Particle acceleration in relativistic collisionless shocks: Fermi process at
  last?", ApJ 682, L5 (2008).
- L. Sironi, A. Spitkovsky and J. Arons, "The maximum energy of accelerated particles in relativistic
  collisionless shocks", ApJ 771, 54 (2013).
- J. G. Kirk, A. W. Guthmann, Y. A. Gallant and A. Achterberg, ApJ 542, 235 (2000).
- A. Achterberg, Y. A. Gallant, J. G. Kirk and A. W. Guthmann, MNRAS 328, 393 (2001).
- J. L. Synge, *The Relativistic Gas* (North-Holland, 1957).
