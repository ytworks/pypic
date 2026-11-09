[English](Magnetic-Reconnection) · **日本語** · ソース: [`examples/reconnection`](https://github.com/ytworks/pypic/tree/main/examples/reconnection)

# 相対論的磁気リコネクション

強く磁化した電子陽電子プラズマ（磁化パラメータ σ = B₀²/(n_b,tot mc²) = {{ c.pair_sigma10.params.sigma }}。
磁場のエネルギー密度 B₀²/2 は背景プラズマの静止質量エネルギー密度の σ/2 = {{ c.pair_sigma10.params.sigma | half }} 倍）の中で、2 枚の薄い電流シートが逆向きの磁場領域を隔てています。シートは磁気島（プラズモイド）の列に
ちぎれ、磁力線がつなぎ変わり、解放された磁場エネルギーが粒子に渡ってべき乗則のスペクトルを作ります。大きな
本番用の計算領域では、2 つの実現で測ったべき指数が p = {{ r.pair_sigma10_large.power_law_index | .2f }} と
{{ r.pair_sigma10_large_run2.power_law_index | .2f }} で、Werner et al.（2016）の σ スケーリングによる
{{ r.pair_sigma10_large.werner2016_index | .2f }} と一致します。

## 物理

**解析的な平衡.** 周期的な二重 Harris シート。半厚さ δ、シートの位置を y₁、y₂ とすると

```math
B_x(y) = B_0\left[\tanh\frac{y-y_1}{\delta} - \tanh\frac{y-y_2}{\delta} - 1\right],
\qquad n(y) = n_{cs}\,\mathrm{sech}^2\frac{y-y_i}{\delta} + n_b .
```

シートの粒子は z 方向に ±β_d でドリフトして電流を担います。アンペールの法則から n_cs β_d = B₀/(2δ)、圧力の
釣り合いから B₀²/2 = 2(n_cs/Γ_d)θ_cs です。σ = B₀²/(n_b,tot mc²) = {{ c.pair_sigma10.params.sigma }} のとき
B₀ = {{ c.pair_sigma10.planner.B0 | .3f }}、β_d = {{ c.pair_sigma10.planner.beta_d | .3f }}、
θ_cs = {{ c.pair_sigma10.planner.theta_cs | .3f }}、アルヴェン速度は v_A = c√(σ/(1+σ)) =
{{ c.pair_sigma10.planner.v_A | .3f }} c です。

**経験的な比較（解析解ではない）.** 相対論的なペアプラズマのリコネクションは約 0.1 v_A の率で進み、粒子
スペクトルは σ に依存する指数のべき乗則 dN/dγ ∝ γ^−p になります。p ≈ 1.9 + 0.7/√σ（Werner et al. 2016）。
Sironi & Spitkovsky（2014）は σ = 1、3、10、50 で p ≈ 4、3、2、1.5 を報告しています。

## 実験条件

規格化単位: c = ε₀ = 1、長さは（背景プラズマの）c/ω_pe、時間は 1/ω_pe。準 2 次元（z 方向 4 セル）、スキン長
あたり 4 セル、全方向周期境界。各粒子種は背景ビーム 1 本とシートのビーム 2 本からなります。補償付き二項電流
フィルタを使います。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![磁気リコネクション](images/reconnection/pair_sigma10_large_run2.gif)

`pair_sigma10_large`（2 つ目の実現、400 × 200 c/ω_pe、t = 2000 まで）。

黒い線は磁力線です。磁束関数 A_z（B = ∇×(A_z ẑ)）の等高線を、全フレームで同じ値の間隔で描いています。
× は磁力線が切れてつなぎ変わる X 点、○ は磁気島（プラズモイド）の中心である O 点です。

- **左上、J_z と磁力線**: 各シートの両側でまっすぐ逆向きに走っていた磁力線が、並んだ X 点で切れ、O 点を囲む
  閉じたループにつなぎ変わります。島はシートに沿って動いて合体し、その間の X 点で磁力線がつなぎ変わり続け、
  最後には各シートが 1〜2 個の大きな島になります。
- **右上、電子密度と磁力線**: プラズマは島に掃き集められ、閉じた磁力線に囲まれた高密度の塊になります。
- **左下、つなぎ変わった磁束 Ψ(t)**（下側のシート）: 主 O 点と最も深い X 点の間の磁束 |A(O) − A(X)| で、磁力線がつなぎ変わるたびに
  増えます。
- **右下、電子のスペクトル**: 初期の熱的分布（灰色）から、γ ≈ 50 以上まで伸びる硬いべき乗則のテイルが育ちます。

## 磁力線とつなぎ変わり

![6 つの時刻の磁力線](images/reconnection/fieldlines.png)

`pair_sigma10_large`（2 つ目の実現）。t = 0、下側シートのつなぎ変わった磁束が最終値の 15、40、65、90 % に
達した時刻、および終了時の J_z と磁力線（全パネルで A_z の同じ等高線値）、下段はそれらの時刻を印した Ψ(t)
です。Ψ = |A(O) − A(X)| は印を付けた点そのもの（下側シートの主 O 点と最も深い X 点）の間で測るので、テアリング
で X 点と O 点ができてから曲線が始まります（それ以前は測る対象がありません）。

- **t = 0**: 磁力線は各シートをはさんでまっすぐ逆向きに走っています。Harris 平衡です。
- **テアリング**: シートが X 点（×）と小さな島（○）の列にちぎれ、X 点で交差していた磁力線は島を囲む閉じた
  ループになります。
- **合体**: 隣り合う島が合体し、その間の磁力線が新しい X 点でつなぎ変わって、Ψ が段階的に増えます（下段）。
- **終盤**: 各シートには少数の大きな島と少しの X 点だけが残り、自由な磁場エネルギーを使い切ると Ψ は頭打ちに
  なります。

## 発表済みの結果との比較

### σ に対するべき指数

![σ に対するべき指数](images/reconnection/power_law.png)

σ = 10 での各 run の測定指数（フィット誤差付きの点。横方向の小さなずれは点を見分けるためだけのもの）、
Werner et al.（2016）のフィット（破線）、Sironi & Spitkovsky（2014）の値（四角）。大きな計算領域は Werner の
曲線の上に乗ります。小さなゲート用領域（100 × 100、t = 400）ではスペクトルが急になります。べき乗則が育つ
には、より大きな系とより長い時間が必要だからです。

### リコネクション率、スペクトル、電流構造

![大きな run](images/reconnection/pair_sigma10_large-reconnection.png)

`pair_sigma10_large`（1 つ目の実現）。**左**: リコネクション率が最大の時刻の J_z。**中央**: example の解析が報告する
リコネクション率 d(Ψ_y₁)/dt/(B₀v_A) と、発表値の ≈ 0.1。Ψ_y₁ = max_x A_z − min_x A_z は、固定した直線 y = y₁
（下側シートの初期位置）に沿って取った値です。磁力線の図の X/O 点間の磁束とは別の指標で、シートが曲がると
固定した直線が膨らんだ磁力線も横切るため、余分に磁束がつなぎ変わらなくても Ψ_y₁ は増えます。**右**: 最終スペクトルと
σ/2 ≤ γ − 1 ≤ 5σ でのべき乗則フィット。

<details><summary>2 つ目の実現とゲート用 run</summary>

![大きな run、2 つ目の実現](images/reconnection/pair_sigma10_large_run2-reconnection.png)
![ゲート用 run](images/reconnection/pair_sigma10-reconnection.png)

</details>

| run | 計算領域 [c/ω_pe] | p（dN/dγ ∝ γ^−p） | Werner 2016 | y = y₁ 上の率の最大値 | y = y₁ 上の活動期の平均率 | 変換された磁場エネルギー | E∥ の寄与（上位 1 %） | 平衡からのずれ / B₀ |
|---|---|---|---|---|---|---|---|---|
| `pair_sigma10`（ゲート） | 100 × 100 | {{ r.pair_sigma10.power_law_index | .2f }} ± {{ r.pair_sigma10.power_law_err | .2f }} | {{ r.pair_sigma10.werner2016_index | .2f }} | {{ r.pair_sigma10.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10.equilibrium_max_dev_over_B0 | .4f }} |
| `pair_sigma10_large` | 400 × 200 | {{ r.pair_sigma10_large.power_law_index | .2f }} ± {{ r.pair_sigma10_large.power_law_err | .2f }} | {{ r.pair_sigma10_large.werner2016_index | .2f }} | {{ r.pair_sigma10_large.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10_large.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10_large.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10_large.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10_large.equilibrium_max_dev_over_B0 | .4f }} |
| `pair_sigma10_large_run2` | 400 × 200 | {{ r.pair_sigma10_large_run2.power_law_index | .2f }} ± {{ r.pair_sigma10_large_run2.power_law_err | .2f }} | {{ r.pair_sigma10_large_run2.werner2016_index | .2f }} | {{ r.pair_sigma10_large_run2.reconnection_rate_peak | .3f }} | {{ r.pair_sigma10_large_run2.reconnection_rate_mean_active | .3f }} | {{ r.pair_sigma10_large_run2.magnetic_energy_converted | .0% }} | {{ r.pair_sigma10_large_run2.E_parallel_share_top1pct | +.3f }} | {{ r.pair_sigma10_large_run2.equilibrium_max_dev_over_B0 | .4f }} |

- **ゲート**（`pair_sigma10`）: テアリングが始まる前（t ≤ 10）、B_x(y) の分布は解析的な平衡から B₀ の 5 % 以内に
  保たれ、ガウスの法則が成り立ちます（{{ r.pair_sigma10.gauss_max | .1e }}）。それ以外はすべて報告のみです。
- **率**: 最大値は初期のテアリングによるものです。活動期の平均は、ゲート用領域では発表値の 0.1 v_A に近く、
  大きな領域では低くなります。大きな領域の率には、初期のピークの後にいくつかのバーストが見えます（上の中央の
  パネル）。example の README は、これをプラズモイドの合体によるものとしています。
- **加速の経路**: 追跡した電子のうち最もエネルギーの高い 1 % は、B に平行な電場成分（E∥）からはほとんど
  エネルギーを得ていません（E∥ の寄与 ≈ 0）。エネルギーは E⊥、つまり外向きの流れとプラズモイドの運動による
  電場から来ています。追跡データは 2.5 ω_pe⁻¹ ごとの記録なので、この内訳は目安です。

## 補足

- **大きな run の 2 つの実現**: 同じ設定・同じシードで、この設定に影響しない 2 つのコミットで実行したもの
  です。GPU の電流割り付けは原子的加算の順序が非決定的で、乱流が丸め誤差の違いを増幅するため、2 つの run は
  独立な実現になります。指数はフィット誤差の範囲で一致し
  （{{ r.pair_sigma10_large.power_law_index | .3f }} と {{ r.pair_sigma10_large_run2.power_law_index | .3f }}）、
  平均率は {{ r.pair_sigma10_large.reconnection_rate_mean_active | .3f }} と
  {{ r.pair_sigma10_large_run2.reconnection_rate_mean_active | .3f }} v_A です。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/reconnection/run.py --config examples/reconnection/configs/pair_sigma10.yaml \
    --out output/wiki/reconnection/pair_sigma10
python examples/reconnection/analyze.py output/wiki/reconnection/pair_sigma10 --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- G. R. Werner, D. A. Uzdensky, B. Cerutti, K. Nalewajko and M. C. Begelman, "The extent of power-law
  energy spectra in collisionless relativistic magnetic reconnection in pair plasmas", ApJL 816, L8 (2016).
- L. Sironi and A. Spitkovsky, "Relativistic reconnection: an efficient source of non-thermal particles",
  ApJL 783, L21 (2014).
- G. R. Werner, D. A. Uzdensky, M. C. Begelman, B. Cerutti and K. Nalewajko, "Non-thermal particle
  acceleration in collisionless relativistic electron–proton reconnection", MNRAS 473, 4840 (2018).
- E. G. Harris, "On a plasma sheath separating regions of oppositely directed magnetic field",
  Nuovo Cimento 23, 115 (1962).
