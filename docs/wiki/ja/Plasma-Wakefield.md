[English](Plasma-Wakefield) · **日本語** · ソース: [`examples/wakefield`](https://github.com/ytworks/pypic/tree/main/examples/wakefield)

# レーザーパルスと電子バンチのプラズマ航跡場

冷たいプラズマを短いドライバーが通り抜けると、その後ろにプラズマ波が残ります。各点は ω_pe で振動し、波の模様は
ドライバーについて進みます。ドライバーは、アンテナから入射するレーザーパルス、または剛体の電子バンチです。箱は
ドライバーを追いかけます（移動ウィンドウ）。弱いドライバー 2 つを線形理論と比べます。軸上の航跡場の振幅は、
パルスの後ろで理論の {{ r.laser_linear.amplitude_ratio | .3f }} 倍、バンチの後ろで
{{ r.beam_linear.amplitude_ratio | .3f }} 倍です。
3 つ目の run では、強いパルス（a₀ = {{ c.laser_bubble.params.a0 }}）が後ろに電子のない空洞を作り、その中に電子が捕捉されます。
空洞の半径は {{ r.laser_bubble.cavity_radius | .2f }} c/ω_pe（Lu et al. 2007 のスケーリングでは
{{ r.laser_bubble.theory.radius | .0f }}）で、捕捉された電子の u_x は最大 {{ r.laser_bubble.trapped.u_max | .0f }} に達します。

## 物理

**線形の航跡場。** イオンが静止した冷たいプラズマでは、電子密度は各点で、ω_pe の強制振動子としてドライバーに
応答します。したがって速度 v のドライバーの後ろの航跡場は、波数 k_w = ω_pe / v の波です。場の単位は
m c ω_pe / e です。

**レーザーパルス。** 場の包絡線が a₀ exp(−t²/τ²) exp(−r²/w²) のパルスは、ポンデロモーティブ力で電子を押します。
その後ろには

```math
E_x = k_w A\, e^{-2r^2/w^2} \cos(t - k_w x) , \qquad
E_r = \frac{4r}{w^2} A\, e^{-2r^2/w^2} \sin(t - k_w x) , \qquad
A = a_0^2\, \frac{\tau}{4} \sqrt{\frac{\pi}{2}}\, e^{-\tau^2/8}
```

が残り、τ = 2 で最大になります（Schroeder et al. 2011 の共鳴するガウス型パルス）。a₀ と w は回折するビームの
その場所での値です。a₀² はプラズマの中では真空中の 1/η 倍になります（η = √(1 − ω_pe²/ω₀²)。パルスは密度の
立ち上がりを通る間、エネルギー流束を保ちます）。比 E_r / E_x = 4 r v / w² は a₀ によりません。

**パルスの速度。** Yee 格子の上では、パルスは格子の群速度

```math
v_g = \frac{\sin(k\,\Delta x)/\Delta x}{\sin(\omega_0 \Delta t)/\Delta t} , \qquad
\frac{\sin^2(\omega_0 \Delta t/2)}{\Delta t^2} = \frac{\sin^2(k\,\Delta x/2)}{\Delta x^2} + \frac{1}{4}
```

で進みます。この設定では {{ c.laser_linear.planner.theory.v_g | .5f }} で、格子がない場合の
√(1 − ω_pe²/ω₀²) = 0.99499 より小さくなります。比較には格子の値を使います。

**電子バンチ。** ピーク密度 n_b、長さ σ_x と σ_r、速度 β の剛体のガウス型バンチの後ろには

```math
E_x = n_b \sqrt{2\pi}\, \sigma_x\, e^{-k_w^2 \sigma_x^2/2}\, R(r) \cos(k_w \xi) , \qquad
R(r) = k_w^2 \int_0^\infty r'\, e^{-r'^2/2\sigma_r^2}\, I_0(k_w r_<)\, K_0(k_w r_>)\, dr'
```

が残ります。ξ = x − βt、k_w = 1/β です（Muggli、arXiv:2007.05226 の式 1 と 3 は β = 1 の場合）。

| 量 | レーザーパルス | 電子バンチ |
|---|---|---|
| ドライバーの速度 | {{ c.laser_linear.planner.speed | .5f }} c | {{ c.beam_linear.planner.speed | .5f }} c |
| 軸上の E_x の振幅 | {{ c.laser_linear.planner.theory.Ex_axis_focus | .3e }}（焦点で） | {{ c.beam_linear.planner.theory.Ex_axis | .4f }} |
| 半径方向の構造 | r = w/2 で E_r / E_x = {{ c.laser_linear.planner.theory.Er_over_Ex_half_waist | .4f }} | R(σ_r) / R(0) = {{ c.beam_linear.planner.theory.R_ratio | .4f }} |
| パルスのレイリー長 | {{ c.laser_linear.planner.theory.rayleigh | .1f }} c/ω_pe | — |

**強いパルス（ブローアウト）。** a₀ ≳ 2 で、ウェストを w = 2√a₀ に合わせたパルスは、後ろのほぼ球形の空洞
（半径 R = 2√a₀）から電子を追い出します。イオンは残ります。軸上では空洞の場は直線で、空洞の中心を x_c として
E_x = (x − x_c)/2 です。したがって最大値は √a₀ です。航跡場はパルスについて速度 1 − 3ω_pe²/(2ω₀²) で進みます。
空洞の後端で捕捉された電子は、位相ずれ長 (2/3)(ω₀/ω_pe)² R の間は加速される側にいて、その間に
(2/3)(ω₀/ω_pe)² a₀ m c² のエネルギーを得ます（Lu et al. 2007）。これらは解ではなくスケーリングなので、強いパルスに
ゲートはありません。

| 量（a₀ = {{ c.laser_bubble.params.a0 }}、ω₀ = {{ c.laser_bubble.params.omega0 }} ω_pe） | Lu et al. 2007 のスケーリング |
|---|---|
| 空洞の半径 | {{ r.laser_bubble.theory.radius | .0f }} c/ω_pe |
| 軸上の E_x の傾き | {{ r.laser_bubble.theory.slope | .1f }} |
| 航跡場の速度 | {{ r.laser_bubble.theory.wake_speed | .3f }} c |
| 位相ずれ長 | {{ r.laser_bubble.theory.dephasing_length | .0f }} c/ω_pe |
| そこまでのエネルギー利得 | {{ r.laser_bubble.theory.energy_gain | .0f }} m c² |

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、プラズマの密度 n₀ = 1、場は m c ω_pe / e。ドライバーは
+x の向きに進みます。

- **プラズマ**: 電子とイオン（質量 1836）を、それぞれセルあたり {{ c.laser_linear.ppc }} 個、同じ位置に
  置きます。そのため、ドライバーが来るまでプラズマに場はありません。セルの中の粒子は 2 × 2 × 2 の格子点に
  置くので、乱数で置いた場合の密度の雑音もありません。密度は傾斜をつけて立ち上げます。
- **レーザーパルス**（`laser_linear`）: a₀ = {{ c.laser_linear.params.a0 }}、τ = {{ c.laser_linear.params.duration }}、
  ウェスト w = {{ c.laser_linear.params.waist }}、ω₀ = {{ c.laser_linear.params.omega0 }} ω_pe、y 方向の直線偏光。
  アンテナは、密度の立ち上がりより手前の真空に置いた電流の面です。焦点は
  x = {{ c.laser_linear.params.focus | .2f }} にあります。レイリー長が経路に比べて長くないので、軸上の a₀² は
  経路に沿って最大 13 % 変わります。理論値は測定する場所で評価します。
- **強いパルス**（`laser_bubble`、記録済みの run）: a₀ = {{ c.laser_bubble.params.a0 }}、τ = {{ c.laser_bubble.params.duration }}、
  w = {{ c.laser_bubble.params.waist }} = 2√a₀、ω₀ = {{ c.laser_bubble.params.omega0 }} ω_pe。焦点は密度の傾斜の上端
  （x = {{ c.laser_bubble.params.focus }}）です。箱の幅は {{ c.laser_bubble.L.1 | .0f }} c/ω_pe、横方向のセルは
  {{ c.laser_bubble.dx.1 }} c/ω_pe です。パルスは 172 c/ω_pe のプラズマを通り、ウィンドウは
  {{ c.laser_bubble.window.window_shift_end }} 回シフトします。
- **電子バンチ**（`beam_linear`）: n_b = {{ c.beam_linear.params.bunch_density }} n₀、
  σ_x = σ_r = {{ c.beam_linear.params.sigma_x }}、γ = {{ c.beam_linear.params.gamma }}。バンチの粒子の質量は
  10⁶ なので減速しません。バンチは、静止した逆電荷のパートナーに重ねて置くので、run は中性で場のない状態から
  始まります。パートナーは後ろに残り、その電荷をプラズマが遮蔽します。
- **移動ウィンドウ**: ドライバーが箱の x 上端の {{ c.laser_linear.params.window_lead }} c/ω_pe 手前に来たときに
  箱が動き始め、その後は光速で 1 セルずつ進みます。x 下端は開いた面です。x 上端は、幅 1 c/ω_pe の真空の余白の
  後ろにある導体壁で、プラズマは箱が届いたところから面ごとに置いていきます。位置と出力は実験室系の値です。
  ウィンドウのシフトは、パルスの run で {{ c.laser_linear.window.window_shift_end }} 回、バンチの run で
  {{ c.beam_linear.window.window_shift_end }} 回です。
- **横方向の面**: パルスでは周期境界、バンチでは導体壁です（周期境界では、パートナーの電荷が航跡場の約 2 % の
  一様な E_x を加えます）。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![レーザーパルスの航跡場](images/wakefield/laser_linear.gif)

軸を通り偏光に垂直な面を、実験室系での位置で示します。枠はウィンドウについて動きます。

- **上、E_x**: 航跡場。
- **中、E_z**: この面での航跡場の半径方向の場。
- **下、E_y**: レーザーパルス。

![電子バンチの航跡場](images/wakefield/beam_linear.gif)

- **上の E_x** と**中の E_z** は上と同じです。
- **下、バンチの密度。**

![強いレーザーパルスの航跡場](images/wakefield/laser_bubble.gif)

強いパルス（`laser_bubble`）を同じ面で示します。

- **1 段目、電子密度**: パルスの後ろの空洞と、それに続く空洞。
- **2 段目、E_x**: E_x が負（青）のところ、つまり空洞の後ろ半分で、電子は前向きに加速されます。
- **3 段目、E_y**: レーザーパルス。
- **4 段目、電子の位相空間 (x, u_x)**: 捕捉された電子。最初の空洞の後端から前に並びます。

## 線形理論との比較

![パルスの航跡場](images/wakefield/laser_linear-wakefield.png)

![バンチの航跡場](images/wakefield/beam_linear-wakefield.png)

最後の出力での、軸を通る面の E_x と、軸上の E_x（当てはめと理論の振幅つき）です。当てはめの式は
E_x = c A_th(x) cos(k x + φ) + b で、A_th は理論の振幅、c、k、φ、b は自由なパラメータです。範囲はドライバーの
後ろで、箱の x 下端から 1 c/ω_pe 以上離れたところです。半径方向の構造は、r = w/2（パルス）または
r = σ_r（バンチ）の 2 本の線の上で、同じ k と φ を使って測ります。

| 量 | PIC | 線形理論 | ゲート |
|---|---|---|---|
| パルス: 軸上の E_x の振幅、PIC / 理論 | {{ r.laser_linear.amplitude_ratio | .4f }} | 1（振幅は焦点で {{ r.laser_linear.predicted_amplitude | .3e }}） | 5 % 以内 |
| パルス: 波数、k v_g | {{ r.laser_linear.k_times_speed | .4f }} | 1 | 2 % 以内 |
| パルス: r = w/2 での E_r / E_x、PIC / 理論 | {{ r.laser_linear.Er_over_Ex_ratio | .4f }} | 1（比は {{ r.laser_linear.Er_over_Ex_theory | .4f }}） | 10 % 以内 |
| バンチ: 軸上の E_x の振幅、PIC / 理論 | {{ r.beam_linear.amplitude_ratio | .4f }} | 1（振幅は {{ r.beam_linear.predicted_amplitude | .4f }}） | 5 % 以内 |
| バンチ: 波数、k β | {{ r.beam_linear.k_times_speed | .4f }} | 1 | 2 % 以内 |
| バンチ: r = σ_r での振幅と軸上の振幅の比 | {{ r.beam_linear.R_ratio | .4f }} | {{ r.beam_linear.R_ratio_theory | .4f }} | 10 % 以内 |
| ガウスの法則、max\|∇·E − ρ\| | {{ r.laser_linear.gauss_max | .0e }}（パルス）、{{ r.beam_linear.gauss_max | .0e }}（バンチ） | — | ≤ 1e-10 |

ゲートの基準は、最初の run の前に決めたものです。

## 強いパルス: 空洞と捕捉された電子

![強いパルスの空洞](images/wakefield/laser_bubble-wakefield.png)

最後の出力（t = {{ r.laser_bubble.t | .1f }}）です。上: 軸を通る面の電子密度に、軸上の E_x（橙）と、最初の空洞の中の
E_x のゼロ点で当てはめた直線（破線）を重ねています。中: 電子の位相空間 (x, u_x)。下: 箱の中の全電子の u_x の
スペクトル。破線より右（u_x > {{ r.laser_bubble.trapped.u_min | .0f }}）の電子を、捕捉された電子として数えます。

| 量 | PIC | Lu et al. 2007 のスケーリング |
|---|---|---|
| 空洞の半径 | {{ r.laser_bubble.cavity_radius | .2f }} | {{ r.laser_bubble.theory.radius | .0f }} |
| 空洞の中のゼロ点での E_x の傾き | {{ r.laser_bubble.Ex_slope | .3f }} | {{ r.laser_bubble.theory.slope | .1f }} |
| パルスの a₀ | {{ r.laser_bubble.a0_measured | .2f }}（出力の中の最大は t = {{ r.laser_bubble.history.t.9 | .0f }} の {{ r.laser_bubble.history.a0.9 | .2f }}） | 焦点で {{ c.laser_bubble.params.a0 }} |
| u_x > {{ r.laser_bubble.trapped.u_min | .0f }} の電子 | {{ r.laser_bubble.trapped.number | .1f }} n₀ (c/ω_pe)³。波長 0.8 µm では {{ r.laser_bubble.trapped_electrons_800nm | .2e }} 個、{{ r.laser_bubble.trapped_charge_pC_800nm | .0f }} pC | — |
| これらの電子の u_x | 平均 {{ r.laser_bubble.trapped.u_mean | .0f }}、最大 {{ r.laser_bubble.trapped.u_max | .0f }} | 位相ずれ長で {{ r.laser_bubble.theory.energy_gain | .0f }} |
| ガウスの法則、run 全体の max\|∇·E − ρ\| | {{ r.laser_bubble.gauss_max | .1e }} | ≤ 1e-10（エンジンの契約） |

- **パルスは形を保ちません。** a₀ は t = {{ r.laser_bubble.history.t.9 | .0f }} で {{ r.laser_bubble.history.a0.9 | .2f }} まで上がり、
  {{ r.laser_bubble.a0_measured | .2f }} に戻ります。空洞の半径は t = {{ r.laser_bubble.history.t.2 | .0f }} で {{ r.laser_bubble.history.cavity_radius.2 | .2f }}、最後の出力で
  {{ r.laser_bubble.cavity_radius | .2f }} です。E_x の傾きは 0.41 から 0.46 の間にあります。
- **電子は t = {{ r.laser_bubble.history.t.3 | .0f }} から run の最後まで捕捉され続ける**ので、スペクトルは広くなります。最大の
  u_x は t = {{ r.laser_bubble.history.t.11 | .0f }} 以降、スケーリングの {{ r.laser_bubble.theory.energy_gain | .0f }} を超えています
  （その時点で {{ r.laser_bubble.history.u_max.11 | .0f }}）。run が終わる時点でも増え続けています。このときパルスは、
  位相ずれ長 {{ r.laser_bubble.theory.dephasing_length | .0f }} c/ω_pe のうち 172 c/ω_pe を進んでいます。
- **軸上の E_x の最大と最小**（{{ r.laser_bubble.Ex_axis_max | .2f }} と {{ r.laser_bubble.Ex_axis_min | .2f }}）は、捕捉された電子の
  位置にある幅 0.06 c/ω_pe のスパイクです。空洞の場ではありません。空洞の場は傾きで測っています。

## 補足

- **ドライバーなし**（この example のテスト run。次の 2 項目とともに
  [`docs/plans/2026-11-11-wakefield.md`](https://github.com/ytworks/pypic/blob/main/docs/plans/2026-11-11-wakefield.md)
  に記録があります）では、どちらの箱も厳密に静止したままです。軸上の E_x は 0 で、平坦部の密度はどのサンプルでも
  n₀ です。格子点に置いた粒子には、ウィンドウが動かす雑音がありません。
- **ウィンドウと静止した箱の比較。** プラズマ全体が入る長さの静止した箱で同じ run をすると、同じ航跡場に
  なります。最後の出力の同じ実験室系の区間で比べると、パルスの run の軸上の E_x の差は最大値の 1.9e-9 以下で、
  バンチの run では保存されたビット（float32 の出力）がすべて一致します。静止した箱のセル数は 2.8 倍（パルス）と
  2.3 倍（バンチ）です。
- **格子点への配置。** セルあたり 8 個の粒子を乱数で置くと、静止した箱ではどちらの振幅もゲートを外れます
  （パルスで −8.3 %、バンチで −6.9 %）。この減少は振動が古いほど大きく、セルあたりの粒子数を増やすと小さく
  なります。原因は特定できていません。ドライバーからさらに離れたところでは、格子点に置いた場合も航跡場は
  減衰します。
- **パルスの run の E_x の曲線の前方にあるさざ波**は、レーザーパルスの位置にあります。当てはめはその後ろから
  始めます。
- **強いパルスの run の ∇·B** は、パルスが入る間は最大 1.6e-11、その後は 1.7e-12 から 2.7e-12 で、エンジンの契約の
  絶対値 1e-12 を超えています。大きさ 50 の場を 0.0196 c/ω_pe のセルで割ったときの丸めの 1 単位は 3.6e-13 です。
  ガウスの法則の残差は run の間に大きくなります（run を 4 等分した各区間の最大は 2.0e-11、2.9e-11、3.6e-11、
  3.7e-11）。軸上の残差の原因は特定できていません。

## 再現方法

```bash
python -m pypic.launch -n 8 examples/wakefield/run.py \
    --config examples/wakefield/configs/laser_linear.yaml --out output/wiki/wakefield/laser_linear
python examples/wakefield/analyze.py output/wiki/wakefield/laser_linear --movie
```

`beam_linear` と `laser_bubble` も同じ方法で実行できます。`--set drive=false` を付けると、ドライバーのない箱を
実行します。

## データの出所

{{ table:provenance }}

## 参考文献

- C. B. Schroeder, C. Benedetti, E. Esarey and W. P. Leemans, "Nonlinear pulse propagation and phase
  velocity of laser-driven plasma waves" (2011), [OSTI 1013067](https://www.osti.gov/servlets/purl/1013067).
- P. Muggli, "Beam-driven systems, plasma wakefield acceleration", CERN Accelerator School 2019,
  arXiv:2007.05226.
- W. Lu, M. Tzoufras, C. Joshi, F. S. Tsung, W. B. Mori, J. Vieira, R. A. Fonseca and L. O. Silva,
  "Generating multi-GeV electron bunches using single stage laser wakefield acceleration in a 3D
  nonlinear regime", Phys. Rev. ST Accel. Beams 10, 061301 (2007).
