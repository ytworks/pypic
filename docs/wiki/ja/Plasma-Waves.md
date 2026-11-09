[English](Plasma-Waves) · **日本語** · ソース: [`examples/plasma_waves`](https://github.com/ytworks/pypic/tree/main/examples/plasma_waves)

# 熱的プラズマ波: 雑音から読む分散関係

静止した一様なマクスウェル分布のプラズマは決して静かではありません。熱揺らぎがプラズマのあらゆる固有
モードを励起するので、場の時空間フーリエ変換 |F(k, ω)|² を取ると分散関係 ω(k) の枝が直接見えます。摂動は
一切与えず、シミュレーションを放置して得たスペクトルを線形理論に重ねます。ゲート対象の枝は、k についての
中央値で {{ r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }}（O モード）から
{{ r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .1% }}（ラングミュア波）の精度で一致します。

## 物理

| 設定 | 磁場 | 比較する枝（理論は `pypic.theory.waves`） |
|---|---|---|
| `unmagnetized` | なし | E_x: ラングミュア波。マクスウェル分布の誘電関数の運動論的な根（kλ_D ≪ 1 ではボーム–グロス ω² = ω_pe² + 3k²v_t²）。E_y, E_z: 光の波 ω² = ω_pe² + c²k²（Yee 格子とリープフロッグ時間積分の数値分散を含む） |
| `perpendicular` | B₀ ∥ z、k ⊥ B₀、ω_ce = {{ c.perpendicular.params.wce }} ω_pe | E_z: 正常波（O モード）。E_y: 異常波（X モード）の上側・下側の枝。E_x: サイクロトロン高調波の間の電子バーンスタイン波 |
| `parallel` | B₀ ∥ x、k ∥ B₀、ω_ce = {{ c.parallel.params.wce }} ω_pe | E_y, E_z: ホイッスラー波、R 波、L 波。E_x: ラングミュア波 |
| `ion_acoustic` | なし、m_i/m_e = {{ c.ion_acoustic.params.mass_ratio }}、T_e/T_i = 50 | E_x: イオン音波。電子とイオンのランダウ減衰を含む運動論的な根 |

冷たい磁化プラズマの枝はスティックスの誘電テンソルから求まります。たとえば X モードは次を満たします。

```math
n^2 = \frac{c^2k^2}{\omega^2} = \frac{(\omega^2-\omega_{pe}^2)^2 - \omega^2\omega_{ce}^2}
{\omega^2(\omega^2 - \omega_{uh}^2)}, \qquad \omega_{uh}^2 = \omega_{pe}^2 + \omega_{ce}^2 ,
```

電子バーンスタイン波は、熱的ラーマー半径 ρ_e = v_te/ω_ce を使った垂直方向の熱い（静電的な）誘電関数の根です。

```math
1 = \frac{2\omega_{pe}^2}{\omega_{ce}^2}\sum_{n=1}^{\infty}
\frac{n^2\, e^{-\lambda} I_n(\lambda)/\lambda}{\omega^2/\omega_{ce}^2 - n^2},
\qquad \lambda = k^2\rho_e^2 .
```

λ → 0 で、第 1 バンドは高域混成周波数 ω_uh に、その他のバンドは高調波 nω_ce に近づきます。

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、θ = T/(m c²)。計算領域は k ∥ x の準 1 次元で
Δx = λ_De、ラインアウトは横方向 8 × 8 セルで平均しています。周波数分解能は Δω = 2π/t_end =
{{ c.unmagnetized.planner.frequency_resolution | .4f }} ω_pe（t_end = 300）です。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![熱雑音からのプラズマ波](images/plasma_waves/unmagnetized.gif)

`unmagnetized`。**上段**: 時間とともに描かれる E_x(x, t) と E_y(x, t)。縦波の E_x はゆっくりしたラングミュア
波の波束が入り混じっています。横波の E_y では、光の波がほぼ c で領域を両方向に横切ります（x–t 面で ±45° の
かすかな筋）。**左下**: その時刻のラインアウト。**右下**: エネルギー。運動エネルギーは一定です。場は
ゼロから始まり、数プラズマ周期のうちに熱揺らぎのレベルまで立ち上がった後、ゆっくり増え続けますが、
運動エネルギーより 2 桁以上小さいままです。

## 線形理論との比較

### 理論曲線を重ねた ω–k スペクトル

色はスペクトルパワーの log₁₀、破線と点線が線形理論です。パワーの尾根はすべて理論曲線の上にあります。

![unmagnetized](images/plasma_waves/unmagnetized-omega_k.png)

`unmagnetized`: E_x のラングミュア波、E_y と E_z の光の波。

![perpendicular](images/plasma_waves/perpendicular-omega_k.png)

`perpendicular`: O モード（E_z）、X モードの上側・下側の枝（E_y）、高調波 nω_ce の間のバーンスタイン
バンド（E_x、kρ_e ≈ 3 まで）。

![parallel](images/plasma_waves/parallel-omega_k.png)

`parallel`: ホイッスラー波、R 波、L 波（E_y, E_z）、ラングミュア波（E_x）。

![ion acoustic](images/plasma_waves/ion_acoustic-omega_k.png)

`ion_acoustic`: 低周波のイオン音波の枝。

### 各枝の誤差

![各枝の誤差の中央値](images/plasma_waves/branch_errors.png)

各 k で、理論周波数の ±15 % の範囲でスペクトルのピークを探します。枝が合格となる条件は、k 点の半分以上で
ピークが見つかり、見つかった点の 90 % 以上が理論から max(3 %, 2Δω/ω) 以内にあることです。塗りつぶしの
マーカーはゲート対象の枝、白抜きは報告のみの枝です。

| 設定 | 枝 | 検出 / k 点 | 誤差の中央値 | 最大誤差 | 許容内の割合 | ゲート |
|---|---|---|---|---|---|---|
| `unmagnetized` | E_x: ラングミュア波（運動論） | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).found }} / {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).points }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ex: Langmuir (kinetic).fraction_within_tolerance | .0% }} | 対象 |
| `unmagnetized` | E_y: 電磁波（Yee） | {{ r.unmagnetized.branches.Ey: EM (Yee).found }} / {{ r.unmagnetized.branches.Ey: EM (Yee).points }} | {{ r.unmagnetized.branches.Ey: EM (Yee).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ey: EM (Yee).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ey: EM (Yee).fraction_within_tolerance | .0% }} | 対象 |
| `unmagnetized` | E_z: 電磁波（Yee） | {{ r.unmagnetized.branches.Ez: EM (Yee).found }} / {{ r.unmagnetized.branches.Ez: EM (Yee).points }} | {{ r.unmagnetized.branches.Ez: EM (Yee).median_rel_error | .2% }} | {{ r.unmagnetized.branches.Ez: EM (Yee).max_rel_error | .2% }} | {{ r.unmagnetized.branches.Ez: EM (Yee).fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_z: O モード | {{ r.perpendicular.branches.Ez: O mode.found }} / {{ r.perpendicular.branches.Ez: O mode.points }} | {{ r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ez: O mode.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ez: O mode.fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_y: X モード上側 | {{ r.perpendicular.branches.Ey: X upper.found }} / {{ r.perpendicular.branches.Ey: X upper.points }} | {{ r.perpendicular.branches.Ey: X upper.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X upper.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X upper.fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_y: X モード下側 | {{ r.perpendicular.branches.Ey: X lower.found }} / {{ r.perpendicular.branches.Ey: X lower.points }} | {{ r.perpendicular.branches.Ey: X lower.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X lower.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ey: X lower.fraction_within_tolerance | .0% }} | 報告のみ |
| `perpendicular` | E_x: バーンスタイン 1 | {{ r.perpendicular.branches.Ex: Bernstein 1.found }} / {{ r.perpendicular.branches.Ex: Bernstein 1.points }} | {{ r.perpendicular.branches.Ex: Bernstein 1.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 1.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 1.fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_x: バーンスタイン 2 | {{ r.perpendicular.branches.Ex: Bernstein 2.found }} / {{ r.perpendicular.branches.Ex: Bernstein 2.points }} | {{ r.perpendicular.branches.Ex: Bernstein 2.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 2.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 2.fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_x: バーンスタイン 3 | {{ r.perpendicular.branches.Ex: Bernstein 3.found }} / {{ r.perpendicular.branches.Ex: Bernstein 3.points }} | {{ r.perpendicular.branches.Ex: Bernstein 3.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 3.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 3.fraction_within_tolerance | .0% }} | 対象 |
| `perpendicular` | E_x: バーンスタイン 4 | {{ r.perpendicular.branches.Ex: Bernstein 4.found }} / {{ r.perpendicular.branches.Ex: Bernstein 4.points }} | {{ r.perpendicular.branches.Ex: Bernstein 4.median_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 4.max_rel_error | .2% }} | {{ r.perpendicular.branches.Ex: Bernstein 4.fraction_within_tolerance | .0% }} | 報告のみ |
| `parallel` | E_y: ホイッスラー波 | {{ r.parallel.branches.Ey: whistler.found }} / {{ r.parallel.branches.Ey: whistler.points }} | {{ r.parallel.branches.Ey: whistler.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: whistler.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: whistler.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_y: R 波 | {{ r.parallel.branches.Ey: R.found }} / {{ r.parallel.branches.Ey: R.points }} | {{ r.parallel.branches.Ey: R.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: R.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: R.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_y: L 波 | {{ r.parallel.branches.Ey: L.found }} / {{ r.parallel.branches.Ey: L.points }} | {{ r.parallel.branches.Ey: L.median_rel_error | .2% }} | {{ r.parallel.branches.Ey: L.max_rel_error | .2% }} | {{ r.parallel.branches.Ey: L.fraction_within_tolerance | .0% }} | — |
| `parallel` | E_x: ラングミュア波（運動論） | {{ r.parallel.branches.Ex: Langmuir (kinetic).found }} / {{ r.parallel.branches.Ex: Langmuir (kinetic).points }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).median_rel_error | .2% }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).max_rel_error | .2% }} | {{ r.parallel.branches.Ex: Langmuir (kinetic).fraction_within_tolerance | .0% }} | — |
| `ion_acoustic` | E_x: イオン音波（運動論） | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).found }} / {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).points }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).median_rel_error | .2% }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).max_rel_error | .2% }} | {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).fraction_within_tolerance | .0% }} | — |

「ゲート」は GPU ゲートの run を指します。ゲート設定は `unmagnetized` と `perpendicular` だけで
（`scripts/gpu_gates.py`）、`parallel` と `ion_acoustic` は報告のみです。`parallel` の E_z の行は E_y と
同じ内容なので、上の図にだけ載せています。

ガウスの法則は {{ r.unmagnetized.gauss_max | .1e }}（`unmagnetized`）と {{ r.perpendicular.gauss_max | .1e }}
（`perpendicular`）、∇·B は {{ r.perpendicular.divB_max | .1e }} の精度で成り立っています。

## 補足

- **X モード下側の枝**（報告のみ）: 高域混成共鳴に近づく枝で、そこでは冷たいプラズマ理論に対する熱的補正が
  小さくなくなるため、表示はしますがゲート対象にはしていません。
- **バーンスタイン第 4 バンド**（報告のみ）: この k 範囲では熱雑音による励起が弱く、ピークが見つかる点が
  半分に届きません。
- **イオン音波**（報告のみ）: この枝は強くランダウ減衰するので、スペクトルのピークが周波数分解能
  （Δω = {{ c.ion_acoustic.planner.frequency_resolution | .4f }} ω_pe）に比べて幅広くなります。誤差の中央値は
  {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).median_rel_error | .1% }} で、3 % の許容内にある点は
  {{ r.ion_acoustic.branches.Ex: ion acoustic (kinetic).fraction_within_tolerance | .0% }}（必要なのは 90 %）
  なので、この設定は解析の基準を満たしていません。
- **磁化ケースで θ_e = 10⁻³ とする理由**: θ_e = 10⁻² では相対論的な質量増加（⟨γ⟩ − 1 ≈ 1.5 θ_e）が
  サイクロトロン高調波をずらし、バーンスタインバンドが非相対論の理論より 1.5〜4 % 低く出ます。
- **出力間隔**: ラインアウトは解析範囲の最高周波数のナイキストレートで書き出し、各区間で時間平均します。
  これが高 k の電磁雑音に対するエイリアシング防止フィルタになります。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/plasma_waves/run.py --config examples/plasma_waves/configs/unmagnetized.yaml \
    --out output/wiki/plasma_waves/unmagnetized
python examples/plasma_waves/analyze.py output/wiki/plasma_waves/unmagnetized --movie
```

## データの出所

{{ table:provenance }}

## 参考文献

- B. D. Fried and S. D. Conte, *The Plasma Dispersion Function* (Academic Press, 1961).
- T. H. Stix, *Waves in Plasmas* (AIP, 1992).
- I. B. Bernstein, "Waves in a plasma in a magnetic field", Phys. Rev. 109, 10 (1958).
