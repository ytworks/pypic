[English](Landau-Damping) · **日本語** · ソース: [`examples/landau_damping`](https://github.com/ytworks/pypic/tree/main/examples/landau_damping)

# ラングミュア振動: ランダウ減衰・Dawson 振動・波の砕波

電子プラズマ（ラングミュア）振動はプラズマのもっとも単純な集団運動で、その性質のいくつかは厳密に分かって
います。この例では 5 つを確かめます。暖かいプラズマ中の小振幅の波の無衝突（ランダウ）減衰、大振幅の波の
非線形減衰と粒子捕捉、Dawson の厳密な冷たい振動、冷たい波が砕ける瞬間、そして相対論的な振動の周期です。
線形ランダウ減衰率は運動論的な根から {{ r.landau.gamma_rel_error | +.1% }} のずれで再現され、3 つの冷たい
厳密解は {{ r.relativistic.freq_rel_error_median | .2% }}（相対論的周期）から
{{ r.wavebreaking.rel_error | abs | .1% }}（砕波時刻）の精度で再現されます。

## 物理

| 設定 | 構成 | 厳密解・理論値 |
|---|---|---|
| `landau` | 暖かい電子（θ_e = {{ c.landau.params.theta_e }}）、n = n₀(1 + α cos kx)、α = {{ c.landau.params.alpha }}、kλ_De = {{ c.landau.params.k_lambda_D }} | マクスウェル分布の誘電関数の運動論的な根: ω = {{ r.landau.omega_theory | .4f }} ω_pe、γ = {{ r.landau.gamma_theory | .4f }} ω_pe |
| `nonlinear` | 同じ構成で α = {{ c.nonlinear.params.alpha }} | 強い減衰のベンチマーク（Cheng & Knorr 1976）: 包絡線はまず γ₁ ≈ −0.29 で減衰し、その後粒子捕捉により γ₂ ≈ +0.08 で再び成長する |
| `dawson` | 冷たい電子、速度 v = v₀ sin kx₀、kv₀ = 0.5 ω_pe | 各電子シートは厳密に x = x₀ + (v₀/ω_pe) sin(kx₀) sin(ω_pe t) で振動する |
| `wavebreaking` | 冷たい電子、kv₀ = 1.2 ω_pe | シートが最初に交差する（波が砕ける）のは ω_pe t = arcsin(ω_pe/kv₀) = {{ r.wavebreaking.t_cross_theory | .3f }} |
| `relativistic` | 冷たい電子、β_max = {{ r.relativistic.beta_max | .2f }} | 不変量 γ + ω_pe²ξ²/(2c²) = 一定 から求まる相対論的なシート周期。T ≈ (2π/ω_pe)(1 + 3β²/16) |

**ランダウ減衰.** 暖かいプラズマでは、波はその位相速度で動く電子とエネルギーをやり取りします。マクスウェル
分布では位相速度より遅い電子のほうが速い電子より多いので、衝突がなくても波はエネルギーを失って減衰します。
複素周波数は次の式の根です。

```math
\varepsilon(k,\omega) = 1 - \frac{1}{2k^2\lambda_{De}^2}\, Z'\!\left(\frac{\omega}{\sqrt{2}\,k v_{te}}\right) = 0 ,
```

ここで Z はプラズマ分散関数（Fried & Conte）です。

**Dawson の冷たい振動.** 冷たいプラズマでは、シートが互いに交差しない限り、振幅によらず各電子シートが ω_pe の
調和振動子になります。kv₀ > ω_pe だと隣のシートが追い越し、密度が無限大になって波は砕けます。

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、θ = T/(m c²)。計算領域は 1 波長の長さで、イオンは
固定された中和背景（m_i/m_e = 1836）です。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

**初期状態.** 密度摂動は、配置した電子を動かして作ります。この移動を電荷保存型のスキームで電流として
割り付けるので、初期の E_x は丸め誤差の範囲でガウスの法則を満たします。陽解法では、その後運動量を Δt/2 だけ
過去に戻します。冷たいケースは、速度摂動だけを与えた場のない平衡状態から始めます。

## アニメーション

### 非線形ランダウ減衰と粒子捕捉（`nonlinear`）

![非線形ランダウ減衰](images/landau_damping/nonlinear.gif)

- **左上、電子の位相空間 (x, u_x)**: 大振幅の波はまず減衰し、位相速度付近の電子が波に捕まって渦（ホール）を
  作り、バウンス周波数で回転します。
- **右上、E_x(x, t)**: 波は最初の ω_pe t ≈ 15 で減衰し、捕捉された電子がエネルギーを返すことで部分的に
  回復します。
- **左下、エネルギースペクトル** と **右下、エネルギー**: 電場のエネルギーが電子に移り、一部が戻ります。

### 波の砕波（`wavebreaking`）

![波の砕波](images/landau_damping/wavebreaking.gif)

追跡した 1,500 個の電子（シート）から描いたもので、色は初期位置 x₀ です。

- **左、位相空間 (x, u_x)**: 冷たい電子は 1 本の曲線上にあります。kv₀ = 1.2 ω_pe では曲線が急になって垂直に
  なり、さらに折り返します。速いシートが遅いシートを追い越すからです。
- **中央、ラグランジュ写像 x₀ → x**: シートの順序が保たれている間、写像は単調です。波が砕ける瞬間に単調で
  なくなります。理論ではこの最初の交差は ω_pe t = {{ r.wavebreaking.t_cross_theory | .3f }} です。
- **右、電子密度**: シートが交差する場所で密度が鋭く立ち上がります（理想的なシートなら無限大）。

交差の後はシートの軌道はもはや Dawson の解ではありません。振動は壊れ、電子は位相空間で混ざり合います。

## 理論との比較

### 許容値に対する誤差

![各比較の誤差](images/landau_damping/errors.png)

各マーカーは 1 つの量の |PIC − 理論|/理論、縦の線はそのゲートの許容値です。白抜きのマーカーは非線形
ベンチマークで、許容値なしの報告のみです。

| 設定 | 量 | PIC | 理論 | 誤差 | ゲート |
|---|---|---|---|---|---|
| `landau` | 減衰率 γ [ω_pe] | {{ r.landau.gamma_fit | .4f }} | {{ r.landau.gamma_theory | .4f }} | {{ r.landau.gamma_rel_error | +.1% }} | ±10 % |
| `landau` | 周波数 ω [ω_pe] | {{ r.landau.omega_fit | .4f }} | {{ r.landau.omega_theory | .4f }} | {{ r.landau.omega_rel_error | +.2% }} | ±2 % |
| `nonlinear` | 初期の減衰 γ₁ | {{ r.nonlinear.gamma_1_fit | .3f }} | −0.29（Cheng & Knorr） | — | 報告のみ |
| `nonlinear` | 捕捉による成長 γ₂ | {{ r.nonlinear.gamma_2_fit | .3f }} | +0.08（Cheng & Knorr） | — | 報告のみ |
| `dawson` | シートの周波数（{{ r.dawson.sheets }} シートの中央値） | — | ω_pe | {{ r.dawson.freq_rel_error_median | .2% }} | 1 % |
| `dawson` | シートの振幅（中央値） | — | v₀ \|sin kx₀\|/ω_pe | {{ r.dawson.amp_rel_error_median | .2% }} | 2 % |
| `wavebreaking` | 最初のシート交差 [1/ω_pe] | {{ r.wavebreaking.t_cross | .3f }} | {{ r.wavebreaking.t_cross_theory | .3f }} | {{ r.wavebreaking.rel_error | +.1% }} | ±5 % |
| `relativistic` | シートの周波数（{{ r.relativistic.sheets }} シートの中央値） | — | 厳密な周期 | {{ r.relativistic.freq_rel_error_median | .3% }} | 0.5 % |

相対論的なケースでは周波数のずれそのものが {{ r.relativistic.nonrelativistic_shift | .1% }}（シートについての
1 − ω/ω_pe の中央値）あるので、測定された周期は相対論的補正を十分な余裕をもって分解しています。ガウスの法則は
全 run で丸め誤差の範囲で成り立っています（max |∇·E − ρ|: `landau` {{ r.landau.gauss_max | .1e }}, `nonlinear` {{ r.nonlinear.gauss_max | .1e }}, `dawson` {{ r.dawson.gauss_max | .1e }}, `wavebreaking` {{ r.wavebreaking.gauss_max | .1e }}, `relativistic` {{ r.relativistic.gauss_max | .1e }}）。

### 線形ランダウ減衰

![減衰するモード（線形）](images/landau_damping/landau-mode.png)

摂動を与えたモードの |E_x(k)| を対数軸で示します。振動は、運動論的な減衰率 γ = {{ r.landau.gamma_theory | .4f }}
ω_pe を傾きとする破線に沿って、熱雑音の床に達するまで減衰します。

### 非線形減衰と粒子捕捉

![減衰するモード（非線形）](images/landau_damping/nonlinear-mode.png)

大振幅の波は、線形理論（破線）より速く、Cheng & Knorr の γ₁（点線）に近い率でまず減衰します。その後、捕捉
された電子によって包絡線が γ₂（一点鎖線）に近い率で再び成長します。

### 冷たいシート: Dawson 振動・砕波・相対論的振動

![Dawson のシート](images/landau_damping/dawson-sheets.png)
![砕波のシート](images/landau_damping/wavebreaking-sheets.png)
![相対論的なシート](images/landau_damping/relativistic-sheets.png)

追跡した電子（シート）の位置 x(t)。`dawson` では各シートが ω_pe の純粋な正弦振動をし、`wavebreaking` では
ω_pe t ≈ 1 の直後に軌道が交差し、`relativistic` では厳密解の予言どおり周期が 2π/ω_pe より長くなります。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/landau_damping/run.py --config examples/landau_damping/configs/landau.yaml \
    --out output/wiki/landau_damping/landau
python examples/landau_damping/analyze.py output/wiki/landau_damping/landau --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- L. D. Landau, "On the vibrations of the electronic plasma", J. Phys. USSR 10, 25 (1946).
- J. M. Dawson, "Nonlinear electron oscillations in a cold plasma", Phys. Rev. 113, 383 (1959).
- C. Z. Cheng and G. Knorr, "The integration of the Vlasov equation in configuration space",
  J. Comput. Phys. 22, 330 (1976).
- B. D. Fried and S. D. Conte, *The Plasma Dispersion Function* (Academic Press, 1961).
