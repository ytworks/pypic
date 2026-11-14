[English](Alfven-Inflow) · **日本語** · ソース: [`examples/alfven_inflow`](https://github.com/ytworks/pypic/tree/main/examples/alfven_inflow)

# 流れに運ばれて衝撃波へ入るアルヴェン的な波

磁化した電子–イオンプラズマ（m_i/m_e = {{ c.gate.params.mass_ratio }}）が、磁場 B₀ = B₀ ŷ を横切って −x
方向へ流れます。プラズマは x+ 面から注入され、B₀ に沿って伝わる円偏光の波を運び込みます。波は冷たい 2 流体
プラズマの R 波で、k d_i = {{ c.gate.planner.k_di | .3f }}、|δB| = {{ c.gate.params.amplitude }} B₀ です。
開いた箱（ゲート）では、波は箱を横切り、その k モードは入射波との差が平均で {{ r.gate.Bx_mean | .1% }}（B_x）、
{{ r.gate.Bz_mean | .1% }}（B_z）です。x− に反射壁を置くと（ピストン）、流れが積み上がって垂直衝撃波
（衝撃波の系でのアルヴェン・マッハ数 {{ r.piston.mach_alfven_shock | .1f }}）ができ、波を運ぶプラズマの中へ
進みます。

この例は、エンジンの 3 つの機能を使います。乱流流入面（PIC-D26）、凍結された初期状態（PIC-D24: ベクトル
ポテンシャル、形を与えた速度、電子の変位による E）、再開用チェックポイント（PIC-D25。ピストンの設定が書き
出します）です。

## 物理

**冷たいプラズマの平行伝搬波。** B₀ = B₀ ŷ の中で静止した電子とイオンは、波 E e^{i(ky − ωt)} に冷たい応答
−iω v = (q/m)(E + v × B₀) で動きます。その電流で、(z, x) 面の円偏光ベクトルについてマクスウェル方程式が
閉じます。サイクロトロン共鳴より下の 2 つの根が、電子と同じ向きに回る R 波と、イオンと同じ向きに回る L 波です。
これらは教科書の関係式と一致します。

```math
\frac{k^2 c^2}{\omega^2} = 1 - \sum_s \frac{\omega_{ps}^2}{\omega(\omega \pm \Omega_s)} ,
\qquad \frac{\omega}{k v_A} \simeq 1 \pm \frac{k d_i}{2} \quad (k d_i \ll 1) .
```

\+ が R 波です。箱の系ではプラズマは −V₀ x̂ で動きます。k は流れに垂直なので、波の位相 k y − ωt は箱の系でも
同じです。電場は E = E_rest − v_d × B となり、運動による電場 −V₀ δB_z ŷ が加わります。

**流入面。** 流入面の表は、波の B、E と各粒子種の速度を、k·r − ωt + φ の cos と sin の成分として持ちます。
入ってくる粒子は、生まれた位置と時刻の波の速度を受け取ります。場は、流入面上の全場・散乱場境界から入ります。

**ピストン衝撃波。** 静止した壁が流れを止め、衝撃波は速度 V_sh で外へ進み、下流は静止します。衝撃波の系での
上流の速度は U₁ = V₀ + V_sh なので、r = n₂/n₁ = U₁/V_sh です。r は M_A = U₁/v_A の垂直 MHD 衝撃波の
圧縮比です（γ = 5/3）。

```math
(2-\gamma) r^2 + \left[\gamma(\beta+1) + (\gamma-1) M_A^2\right] r - (\gamma+1) M_A^2 = 0 ,
\qquad V_{sh} = \frac{V_0}{r - 1} .
```

下流の温度 T_e + T_i は運動量流束から決まります。プラズマに凍結した磁場なら、法線成分 δB_x は変わらず、
接線成分 δB_z は r 倍に圧縮されます。

| 量 | ゲート | ピストン |
|---|---|---|
| アルヴェン速度 v_A | {{ c.gate.planner.v_a | .4f }} c | {{ c.piston.planner.v_a | .4f }} c |
| 流れの速度 V₀ | {{ c.gate.planner.v0 | .4f }} c | {{ c.piston.planner.v0 | .4f }} c |
| 波の振動数 ω / (k v_A)（R 波） | {{ c.gate.planner.omega_over_kva | .3f }} | {{ c.piston.planner.omega_over_kva | .3f }} |
| プラズマ β | {{ c.gate.planner.beta | .2f }} | {{ c.piston.planner.beta | .2f }} |
| ジャンプ条件: r、V_sh、M_A | — | {{ r.piston.compression_theory | .2f }}、{{ r.piston.shock_speed_theory | .4f }} c、{{ r.piston.mach_alfven_shock | .2f }} |

## 実験条件

規格化単位: c = ε₀ = 1、長さは d_e = c/ω_pe、時間は 1/ω_pe、各粒子種の密度 n₀ = 1。z 方向は 2 セルの
準 2 次元です。

- **ゲート（`gate`）**: x− はプラズマが出ていく開放面です。t = 0 で箱全体に波を置きます（B_x、B_z、E_x は
  y 方向のプロファイル、E_z は電子の z 方向の変位、速度は形を与えた摂動）。x+ 面はマーダー領域つきの `flux`
  面です。
- **ピストン（`piston`）**: x− は反射壁です。波に包絡線 h(x) を掛けて置きます。h は壁から 4 d_e 以内で 0、
  12 d_e で 1 になり、流入面の手前 9 d_e で 0 に戻ります。B はベクトルポテンシャル A_z = h a_x(y)、
  A_x = −h a_z(y) から作り（∇·B = 0 が厳密に成り立ちます）、E_x と E_z は電子の変位で作ります。流入面は
  t = 0 から波を運び込みます。4000 ステップごとに再開用チェックポイントを書きます。
- **ヌル run（`piston_null`）**: 波のないピストン（振幅 0 の表）。衝撃波自身の k モードを測ります。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

ゲート: 開いた箱を流れに運ばれて横切る波。

![ゲート](images/alfven_inflow/gate.gif)

ピストン: 衝撃波が波の中へ進みます。

![ピストン](images/alfven_inflow/piston.gif)

- **左、z = 0 の面の B_z / B₀。**
- **中央、z = 0 の面のイオン密度。**
- **右、y と z で平均した電子とイオンの密度。**

## ゲート: 入射波との比較

![ゲート](images/alfven_inflow/gate-alfven_inflow.png)

終了時の B_z、各 x 列の k モードと入射波の比較、密度です。

| 量 | PIC | ゲート |
|---|---|---|
| B_x の列ごとの誤差の最大値（x > 10 d_e） | {{ r.gate.Bx_row_max | .1% }} | ≤ 12 % |
| B_x の列ごとの誤差の平均（10 < x < 60 d_e） | {{ r.gate.Bx_mean | .1% }} | ≤ 4 % |
| B_z の列ごとの誤差の最大値 | {{ r.gate.Bz_row_max | .1% }} | ≤ 12 % |
| B_z の列ごとの誤差の平均 | {{ r.gate.Bz_mean | .1% }} | ≤ 4 % |
| 2 d_e ごとの密度 | {{ r.gate.density_deviation | .1% }} 以内 | 10 % |
| ガウスの法則、max\|∇·E − ρ\| | {{ r.gate.gauss_max | .0e }} | ≤ 1e-10 |
| max\|∇·B\| | {{ r.gate.divB_max | .0e }} | ≤ 1e-12 |
| 平均振幅比 B_x、B_z | {{ r.gate.Bx_ratio | .3f }}、{{ r.gate.Bz_ratio | .3f }} | 報告のみ |

基準は最初の run の前に決めました。列ごとの誤差が最大になるのは流入面のすぐ内側です。

## ピストン: 衝撃波と波

![ピストン](images/alfven_inflow/piston-alfven_piston.png)

終了時の B_z と衝撃波の位置（破線）、x–t 面のイオン密度と衝撃波の位置、衝撃波の前・遠い上流・後ろでの B_x と
B_z の k モード、電子のスペクトルです。これらは報告のみで、ゲートの基準ではありません。

| 量 | 波あり | ヌル run（波なし） | ジャンプ条件 |
|---|---|---|---|
| 衝撃波の速度（後半で当てはめ） | {{ r.piston.shock_speed | .4f }} c | {{ r.piston_null.shock_speed | .4f }} c | {{ r.piston.shock_speed_theory | .4f }} c |
| 下流の密度 n₂ / n₁ | {{ r.piston.compression | .2f }} | {{ r.piston_null.compression | .2f }} | {{ r.piston.compression_theory | .2f }} |
| 下流の B_y / B₀ | {{ r.piston.field_compression | .2f }} | {{ r.piston_null.field_compression | .2f }} | {{ r.piston.compression_theory | .2f }} |
| 下流の (T_e + T_i) / ジャンプ条件 | {{ r.piston.temperature_ratio | .2f }} | {{ r.piston_null.temperature_ratio | .2f }} | 1 |
| 下流の T_e、T_i [m_e c²] | {{ r.piston.Te_down | .2e }}、{{ r.piston.Ti_down | .2e }} | {{ r.piston_null.Te_down | .2e }}、{{ r.piston_null.Ti_down | .2e }} | 和 {{ r.piston.temperature_theory | .2e }} |
| B_x モード 前 / 後ろ（波の振幅単位） | {{ r.piston.Bx_ahead_mean | .2f }} / {{ r.piston.Bx_behind_mean | .2f }} | {{ r.piston_null.Bx_ahead_mean | .3f }} / {{ r.piston_null.Bx_behind_mean | .2f }} | 凍結: 後ろ = 前 |
| B_z モード 前 / 後ろ（波の振幅単位） | {{ r.piston.Bz_ahead_mean | .2f }} / {{ r.piston.Bz_behind_mean | .2f }} | {{ r.piston_null.Bz_ahead_mean | .3f }} / {{ r.piston_null.Bz_behind_mean | .2f }} | 凍結: 後ろ = r × 前 |
| ガウスの法則、max\|∇·B\| | {{ r.piston.gauss_max | .0e }}、{{ r.piston.divB_max | .0e }} | {{ r.piston_null.gauss_max | .0e }}、{{ r.piston_null.divB_max | .0e }} | ≤ 1e-10、≤ 1e-12 |

モードの値は、両方の窓に列がある出力についての平均です。前は x_sh + 4 から x_sh + 20 d_e（反射イオンの
フットがある範囲）、後ろは 0.2 x_sh から x_sh − 4 d_e です。

## 補足

- **衝撃波の後ろの k モードは、凍結した磁場の値より大きくなります。** 衝撃波を横切って B_x は
  {{ r.piston.Bx_transmission | .1f }} 倍、B_z は {{ r.piston.Bz_transmission | .1f }} 倍になります（凍結なら 1 倍と
  {{ r.piston.compression_theory | .2f }} 倍）。ヌル run では、衝撃波自身のその k の揺らぎが、波の振幅単位で
  平均で {{ r.piston_null.Bx_behind_mean | .2f }}（B_x）と {{ r.piston_null.Bz_behind_mean | .2f }}（B_z）、終了時に
  {{ r.piston_null.Bx_behind_end | .2f }} と {{ r.piston_null.Bz_behind_end | .2f }} です（波ありでは
  {{ r.piston.Bx_behind_end | .2f }} と {{ r.piston.Bz_behind_end | .2f }}）。どちらも
  1 回の run で、この例では衝撃波自身のさざ波と、衝撃波を通った波とを分けられません。
- **衝撃波の速度と圧縮比は 2 つの run で異なります**（ジャンプ条件に対して波ありで
  {{ r.piston.shock_speed_ratio | .3f }} 倍、波なしで {{ r.piston_null.shock_speed_ratio | .3f }} 倍）。どちらも
  1 回の run を 1 つの時刻で測った値です。
- **包絡線を流入面でも 0 にする理由。** エンジンは、非周期面に接するベクトルポテンシャルがその面で 0 でないと
  受け付けません（反射壁は法線方向の B を 0 にするため）。そのため流入面の手前 9 d_e は波なしで始まります。
  このプラズマは V₀ で押し流され、t ≈ 850 ごろ衝撃波に届き、その後は流入面から入った波が衝撃波に当たります。
- **ヌル run が表を持つ理由。** エンジンが壁に接する運動による電場 E₀ を許すのは、同じ軸の開放面が揺らぎの表を
  持つときだけです。ヌル run は振幅 0 の表を持ちます。
- 波の y 方向の運動による電場（−V₀ δB_z）は置いていません。その電荷 v₀ δv / c² を粒子が持たないためで、
  プラズマが自分でこの電場を作ります。

## 再現方法

```bash
python -m pypic.launch -n 2 examples/alfven_inflow/run.py \
    --config examples/alfven_inflow/configs/gate.yaml --out output/wiki/alfven_inflow/gate
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml
python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/piston.yaml \
    --set wave=false --out output/alfven_inflow/piston_null
python examples/alfven_inflow/analyze.py output/alfven_inflow/piston --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- T. H. Stix, *Waves in Plasmas* (AIP, New York, 1992), chapter 2: 冷たいプラズマの R 波と L 波。
- E. R. Priest, *Solar Magnetohydrodynamics* (Reidel, Dordrecht, 1982): 垂直 MHD 衝撃波のジャンプ条件。
