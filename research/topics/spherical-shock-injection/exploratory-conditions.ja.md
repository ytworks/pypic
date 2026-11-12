# 試行計算：エンジンで可能な範囲で超新星残骸に近づけた 3 条件

テーマ `spherical-shock-injection`。2026-11-12 設計。**状態：設計し、planner で確認済み。実行はしていない。**
これは `exploratory-conditions.md`（英語、原本）の日本語版である。食い違いがあれば英語版に従う。

**2026-11-13 訂正：** この条件のレビュー（イシュー #48、`reviews/2026-11-13-exploratory-conditions.md` に保存）を
受けて、2.2 節と 2.3 節を直した。残骸の各値に、仮定した速度と電子温度を書いた。「合わせられたもの」は L の β
だけである。L、M、H の物理条件は変えていない。

citekey は `research/literature/` のノートを指す。「調査」は `survey-snr-electron-injection.md` である。
「本記録の計算」はここで導いた値、「見積もり」はどの計算でも確かめていない外挿である。

## 1. 目的と状態

- **目的。** これらの計算は、球面衝撃波で電子がどう加速されるかについて仮説を立てるための試行である。
  主張の証拠ではない。解析のあとで、テーマの研究の問いを絞る（オーナーの判断）。
- **オーナーの依頼（2026-11-12）：**
  - 条件はできるだけ超新星残骸（SNR）に近づける。
  - 空間は大きめ、粒子数は少なめ、計算時間は短く。
  - 2〜3 本。
  - この段階では実行しない。
- **オーナーの決定（2026-11-12）：**
  - 既存の例 `examples/blast_wave` で表せる範囲に限る。新しいコードは書かない。
  - 3 本の違いは磁場の強さだけにする。
  - 格子は 768³、粒子は 1 セル・1 粒子種あたり 2 個。
- **いまあるもの：** 例の planner が受理する設定ファイル 3 本（`configs/`、10 節）と、この記録。GPU でも
  それ以外でも、計算は 1 本も実行していない。
- **最も重要な制約：** この種の爆風の衝撃波面は破片より遅く、しかも減速する。H は名目上サーフィンの閾値より
  上にあるが、面そのものは閾値と同じか、それより下になる見込みである（4.3 節）。

## 2. 超新星残骸にどこまで近いか

### 2.1 何が選択を縛るか
エンジンの単位（c = 1、長さは d_e = c/ω_pe）と例のパラメータでは、次の関係がある。

- M_A = (V_d/c) · (ω_pe/Ω_ce) · (m_i/m_e)^(1/2)。T_i = T_e のとき β = 4 (v_te/c)² (ω_pe/Ω_ce)²
  （`examples/blast_wave/plan.py`）。
- **デバイ長。** 格子幅は 0.5 d_e で、エンジンの実績ではデバイ長を 0.2 格子、つまり v_te = 0.1 c にしている。
  電子をこれより冷たくすると数値加熱が起きる（4.3 節）。v_te を固定すると、同じ流速で M_A を上げることは
  磁場を弱めることであり、β が上がる。
- **箱。** planner は半幅が 1.3 R_B 以上であることを求める。R_B = R_m M_A^(2/3) は磁場による停止半径である。
  破片が磁場を通して周囲のプラズマと結合するのは、破片のジャイロ半径 M_A d_i が等質量半径 R_m を超えない
  場合に限られる（Le et al. 2021。`examples/blast_wave/README.md` の引用による）。768³ の格子と
  m_i/m_e = 16 では、これで M_A は約 8.7 までになる。
- **計算の長さ。** 1 本の長さは 1.25 R_B/V_d で、イオンのジャイロ時間では 1.25 (R_m/d_i) M_A^(−1/3) になる。
  M_A が高いほど、ジャイロ時間で測った長さは短い。

### 2.2 残骸との数値の比較

計算のマッハ数は名目値である：衝撃波面の速度ではなく、破片のドリフト速度 V_d を使っている。

| 量 | 超新星残骸 | L | M | H |
|---|---|---|---|---|
| アルフベン・マッハ数 M_A | 増幅前の磁場で約 250–300、増幅後の磁場で約 15（5000 km/s） | 3.5 | 6.0 | 8.7 |
| M_A / (m_i/m_e)^(1/2)。ホイッスラー波による注入は ≲ 1 が必要 | 増幅前 6–7、増幅後 0.35 | 0.88 | 1.5 | 2.2 |
| M_A / [1.2 (m_i/m_e)^(2/3)]。サーフィンは ≥ 1 が必要 | 増幅前 1.4–1.7、増幅後 0.08 | 0.46 | 0.79 | 1.14 |
| (M_A/√β_e) / [(m_i/m_e)^(1/2)/2] = √2 ×（流速 / 電子の熱速度）。反射電子によるホイッスラー波は ≳ cos θ_Bn が必要 | 5000 km/s・10⁴ K で約 18、3000 km/s・T_e = 10 eV で 3.2 | 3.5 | 3.5 | 3.5 |
| 音速マッハ数 | 約 440（Tycho、4990 km/s） | 5.5 | 5.5 | 5.5 |
| 高速磁気音波マッハ数 | — | 2.95 | 4.05 | 4.64 |
| β | 10⁴ K で 0.4–0.8、10 eV で約 5–10 | 0.49 | 1.44 | 3.03 |
| 衝撃波速度 / c | 0.007–0.02 | 0.25 | 0.25 | 0.25 |
| 球面のうち超光速になる部分 | 2 % 未満 | 25 %（θ_Bn > 76°） | 25 % | 25 % |
| 流速 / 電子の熱速度 | 5000 km/s・10⁴ K で約 13、3000 km/s・T_e = 10 eV で 2.3 | 2.5 | 2.5 | 2.5 |
| m_i/m_e | 1836 | 16 | 16 | 16 |
| ω_pe/Ω_ce | 約 350 | 3.5 | 6.0 | 8.7 |
| 計算の長さ（イオンのジャイロ時間） | — | 13.2 | 7.7 | 5.3 |
| 破片のジャイロ半径 / R_m | — | 0.22 | 0.54 | 1.00 |

残骸の列の出典：
- M_A：`morlino2012strong` Sect. 2.3、`riquelme2011electron` Sect. 5、調査 2 節。約 5000 km/s での値である。
- 閾値：`riquelme2011electron` Sect. 3.1（ホイッスラー波）、`matsumoto2012electron` Eq. 8（サーフィン）、
  `amano2010critical` Eqs. 2, 3（反射電子によるホイッスラー波）。Eq. 3 は
  M_A/√β_e ≃ 68 (V_s / 3000 km/s) (T_e / 10 eV)^(−1/2) で、実際の質量比での閾値は 21.4 である。したがって
  この行は、論文の 3000 km/s・10 eV では 3.2、β の行と流速の行が仮定する 5000 km/s・10⁴ K（0.86 eV）では
  約 18 になる（本記録の計算）。
- この行と「流速 / 電子の熱速度」の行は同じ量である。β_e = 2 (m_e/m_i) (v_te/v_A)²、v_te = (kT_e/m_e)^(1/2)
  なので M_A/√β_e = (m_i/m_e)^(1/2) V/(√2 v_te) となり、この行はどの質量比でもその比の √2 倍である
  （本記録の計算）。
- 音速マッハ数：`morlino2012strong` Sect. 2.3（4990 km/s で M_0 = 443）。
- β：イオンと電子がともに 10⁴ K の場合の本記録の計算。Tycho で 0.8（調査 2 節）、SN 1006 の密度と磁場で 0.4。
  10 eV では 11.6 倍になる。
- ω_pe/Ω_ce：調査 2 節の計算。
- 衝撃波速度：`vink2012supernova` Sect. 11.3。超光速の部分：調査 2 節。
- 流速と電子の熱速度の比：本記録の計算。v_te は 10⁴ K で 390 km/s、10 eV で 1330 km/s である。

L、M、H の値は、planner の出力（10 節）と、そこからの本記録の計算である。

### 2.3 合わせられたもの、合わせられなかったもの
- **合わせられたもの：** L の β。上流温度が 10⁴ K の残骸に対してである。残骸の数値のうち、合っているのは
  これだけである。
- **種類として似ているが、数値を合わせたわけではないもの：**
  - 初期の幾何：一様磁場の中の球対称な駆動なので、動径方向と磁場のなす角は 0° から 90° まで変わる。
    これは設定の説明であり、衝撃波の説明ではない。極の近くに衝撃波ができる見込みはなく（4.2 節）、θ_Bn は
    衝撃波がある場所で、測った面から決めなければならない。
  - ホイッスラー波の上限とサーフィンの閾値に対して、名目の M_A がどちら側にあるか。L はホイッスラー波の
    上限より下にあり、上流磁場が増幅された残骸と同じ側にある。H は名目上サーフィンの閾値より上にあり、
    増幅されていない残骸と同じ側にある。M はその間にある。面そのものは名目値より遅い（4.3 節）。
- **どの計算でも合わせられないもの：**
  - マッハ数そのもの。とくに音速マッハ数。
  - 質量比と衝撃波速度。超光速の帯が 10 倍以上広い。
  - 流速と電子の熱速度の比。`amano2010critical` の閾値に対する M_A/√β_e も同じ量である（2.2 節）。計算は
    2.5（閾値を単位にすると 3.5）。5000 km/s・10⁴ K の残骸は約 13（18）で、β が 0.4–0.8 になるのはこの温度で
    ある。2.3（3.2）になるのは 3000 km/s・T_e = 10 eV の場合だけで、そのとき β は約 10 倍になる。流速と温度を
    固定しているので、この値は 3 本とも同じである。
  - 衝撃波の厚みを単位とした系の大きさ。
- **M と H で合わせられないもの：** β。温度を固定しているので 1.4 と 3.0 まで上がる。したがって、M_A の走査は
  β の走査でもある。

## 3. 3 本の条件

### 3.1 3 本に共通する部分
- **駆動：** `examples/blast_wave` の破片の殻。0.33 R_m にあるガウス型の殻で、幅 0.1 R_m、ピーク密度 16.5 n₀、
  動径方向のドリフト V_d = 0.25 c。粒子種ごとに 4 本のビームとして入れる。
- **周囲のプラズマ：** 静止、T_i = T_e、v_te = 0.1 c、m_i/m_e = 16。6 面の開放境界から注入する。
- **磁場：** z 方向の一様な B₀（外部磁場）。吸収する球体と双極子は置かない（`body: false`）。残骸には磁化した
  天体がないためである。
- **格子：** 0.5 d_e の格子が 768³。一辺は 384 d_e = 96 d_i。dt = 0.25/ω_pe、2960 ステップ。
- **粒子：** 1 セル・1 粒子種あたり 2 個。周囲のプラズマが 1.81 × 10⁹ 個、破片が 1.2 × 10⁷ から 7.5 × 10⁷ 個。
- **数値：** 陽解法、`precision: compact`、補償つきのフィルタ 4 回、Marder ゾーン 0.125、乱数の種 12345。
  8 GPU でのランクの分割は (1, 2, 4)。
- **R_B** は 3 本とも約 148 d_e = 37 d_i（`rm` で調整）。そのため、箱と、1/ω_pe で測った計算の長さは同じになる。

### 3.2 3 本の違い

| 計算 | 設定ファイル | `mach_alfven` = `wpe_wce` | `rm` (d_i) | R_m (d_e) | v_A / c | β | 計算の長さ (1/Ω_ci) |
|---|---|---|---|---|---|---|---|
| L | `configs/explore-L.yaml` | 3.5 | 16.0 | 64 | 0.0714 | 0.49 | 13.2 |
| M | `configs/explore-M.yaml` | 6.0 | 11.2 | 44.8 | 0.0417 | 1.44 | 7.7 |
| H | `configs/explore-H.yaml` | 8.7 | 8.7 | 34.8 | 0.0287 | 3.03 | 5.3 |

### 3.3 なぜ磁場の強さを変えるのか
- ドリフト速度と温度が 3 本とも同じなので、音速マッハ数は同じで、変わるのは磁場だけである。これは、上流磁場が
  増幅された残骸の衝撃波と、増幅されていない衝撃波の違いにあたる（調査 4.1 節：ホイッスラー波は M_A ≲ 20、
  サーフィンは捕捉の閾値より大きい M_A が必要）。
- 調査は、SNR 衝撃波の注入にどのマッハ数が必要かを、文献の食い違いとして挙げている（8.1 節の 2）。3 本は、
  名目のマッハ数で見て、二つの閾値の両側とその間に置いてある。
- M_A が高いほど電子を冷たくすれば、β を揃えられた。それを選ばなかったのは、デバイ長が 0.1 格子まで下がり、
  上流の温度を保てなくなるからである（4.3 節）。

## 4. この計算で見えること、見えないこと

### 4.1 見えること（仮説のための最初の観察）
- 球のまわりの衝撃波面と磁場の形。測った面の法線から求めた局所的な θ_Bn。
- 高エネルギーの電子が面のどこに現れるか。θ_Bn と時間への依存と、L から H への変化。
- 子午面と赤道面での、場所ごとの電子スペクトル。
- 追跡した電子に沿った、エネルギーの履歴と、磁場に平行・垂直な電場の仕事。
- 電子が面より前に出るかどうか、出るならどの θ_Bn か。

### 4.2 見えないこと
- **Buneman 波でのサーフィン。** 波長あたり約 6 格子しかなく（2π · 2V_d/ω_pe ≈ 3.1 d_e、本記録の計算）、
  フィルタで一部が減衰し、1 セル 2 個の雑音に埋もれる。そのような波が現れるかどうかは見えるが、電子を
  どれだけ捕捉するかは見えない。この計算でサーフィンが見えなくても、何も言えない（調査 7 節）。
- **準平行衝撃波。** 準平行衝撃波でイオンが駆動する波には数十ジャイロ時間、フォアショックには約 125 以上が
  必要である（`morris2023preacceleration` Sect. 4.5）。計算の長さは 13、7.7、5.3 である。B₀ に沿っては
  破片が自由に流れ（`examples/blast_wave/README.md`）、計算の終わり近くで z の面に届く。極の近くで見えるのは、
  衝撃波ではなく、プラズマに突っ込むビームである。
- **拡散衝撃波加速。** 箱が小さすぎる（閾値は調査 2 節）。
- **M_A を固定したときの、質量比・衝撃波速度・β への依存。** 変えているのは磁場の強さだけで、β もそれと
  一緒に変わる。
- **θ_Bn > 76° での電子の反射。** 0.25 c ではこの帯は超光速である。面が減速すると帯は狭くなる。

### 4.3 リスク
- **面は V_d より遅く、減速する。** 例の gate 計算では、面は 0.3 から 0.7 R_B の間を 0.87 V_d で動き
  （`docs/plans/2026-11-11-blast-wave.md`、gate の結果）、その先ではもっと遅かった。
  - 0.87 V_d なら、H の面は M_A = 7.6 で、サーフィンの閾値と同じであり、上ではない（本記録の計算）。
  - 面は R_m のスケールで減速する。R_B は H で 4.2 R_m、その gate 計算で 2.1 R_m である。したがって、H の面は
    計算の大部分で閾値をはっきり下回る見込みである。外挿すると、面の M_A は L で約 3、M で 4–5、H で
    高々約 5 になる（見積もり）。
  - 各計算は、名目値ではなく、測った面のマッハ数で呼ぶことにする。
- **H の結合。** 破片のジャイロ半径が R_m と等しい。Le et al. 2021 は、結合の条件をジャイロ半径 < R_m と
  している。例の試行で 0.75 R_m にしたときは、破片のイオンが面の後ろにとどまらなかった
  （`docs/plans/2026-11-11-blast-wave.md`、ビーム 1 本と 4 本の試行）。H では、約 2 R_m = 70 d_e より内側の
  高エネルギー電子は、衝撃波ではなく、周囲のプラズマを突き抜ける破片から来ている可能性がある。
- **上流電子の数値加熱。** 例の null 計算（1 セル 8 個）は、1760 ステップのあと、T_e が内部で 1.5 %、面の
  近くの最悪の板で 3.0 % 上がり、全エネルギーが 2.2 % 増えて終わった（同じ計画、null 試験 N1）。加熱が
  粒子数に反比例するなら、今回の計算は、T_e が内部で約 10 %、面の近くで 13–20 % 上がって終わる（見積もり）。
  これは各計算で測る（6 節）。
  - その場合、M と H では、周囲のプラズマが数値的に得るエネルギーが、破片のエネルギーと同程度か、それより
    大きくなる（見積もり）。例のエネルギー収支の指標は破片のエネルギーを単位にしているので、そのために
    大きな値に見える。
- **雑音。** 1 セル 2 個では場の雑音が大きい。スペクトルの尾には粒子数が要る（6 節）。
- **計算の最後の 1 割。** gate 計算の面が最後に達した 1.17 R_B は、境界の面から 19 d_e の位置である
  （本記録の計算）。各計算の終わりは境界の影響を受ける。
- **試されていない出力の経路。** エネルギーによる粒子の選別（`select`）は、リポジトリにテストも使用例もない。
  また、ステップ 0 のように何も選ばれないダンプは、空のデータセットを書く。最初の計算の前に、GPU の計算機で
  16 ステップの確認をする（8 節）。

## 5. 出力

planner の既定の出力には粒子データがない。そのため、各設定ファイルで診断の一覧を指定している。最初の三つは
planner の診断を同じ名前で繰り返したものである。`examples/blast_wave/analyze.py` がそれを読むためである。

| 名前 | 種類 | 内容 | 間隔（ステップ） | 回数 | サイズ（圧縮前） |
|---|---|---|---|---|---|
| `coarse` | fields | 4 粒子種のモーメント、E、B、ρ。13 / 9 / 7 格子のブロック（L / M / H） | 296 | 11 | 0.26 / 0.76 / 1.6 GB |
| `slice_z` | slice | 赤道面（全解像度）：密度、イオンと電子の流束、電子のエネルギー密度、E、B、J | 74 | 41 | 1.9 GB |
| `slice_y` | slice | 子午面（B₀ を含む面、全解像度）、同じ場 | 74 | 41 | 1.9 GB |
| `spec_meridional` | histogram | 子午面を挟む厚さ 16 d_e の板の中の 4 粒子種：(x, z, 運動エネルギー)、192 × 192 × 96 ビン、エネルギーは 10⁻⁴ から 10² m_e c²、1 桁あたり 16 ビン | 148 | 21 | 1.2 GB |
| `spec_equatorial` | histogram | 赤道面を挟む板で同じもの：(x, y, 運動エネルギー) | 148 | 21 | 1.2 GB |
| `spec_all` | histogram | 箱全体のエネルギースペクトル、4 粒子種 | 74 | 41 | 1 MB 未満 |
| `spec_up` | histogram | 爆風が届かない隅の柱の中の、周囲の粒子種のエネルギースペクトル | 74 | 41 | 1 MB 未満 |
| `hot_e` | particles | γ ≥ 1.5 のすべての電子（周囲と破片）：位置、運動量、重み、id | 148 | 21 | 下記 |
| `tracks_e` | tracks | 周囲の電子 1.9 × 10⁵ 個（id が 4733 で割り切れるもの）を最初から追跡。粒子位置の E と B つき | 4 | 741 | 6.8 GB |

- **`hot_e`：** γ = 1.5 は運動エネルギー 0.5 m_e c² で、イオンのラム・エネルギー（½ m_i V_d²）にあたり、
  上流の電子温度の 50 倍である。
  - 静止した電子が 1 回ミラー反射されると、θ_Bn = 0°、45°、60°、70° で 0.13、0.29、0.67、2.3 m_e c² を得る
    （本記録の計算、相対論的）。したがって γ ≥ 1.5 は、56° 以上での 1 回の反射か、それより小さい角度での
    複数回の相互作用を意味する。
  - planner の表では 1 回あたり 16 GB と出る。これは「全粒子の 1 割」という planner の一律の見積もりである。
    実際の数は、下流の電子温度によって、1 回あたり 10² から数 10⁶ 個の見込みである（見積もり）。
  - 選ばれた粒子はすべて 1 台の GPU に集められる。1 回のダンプが約 10⁷ 個を超える場合は、閾値を γ ≥ 2 に
    上げる。
- **`tracks_e`：** 追跡する集合は最初の出力で固定される。
  - 間引きの 4733 は素数で、追跡する電子が格子の上で規則的に並ばない。768 と約数を共有する数にすると、
    電子が平面の上に並んでしまう。
  - 追跡する電子のうち約 7 × 10⁴ 個が、面の通る範囲にある（本記録の計算）。千個に 1 個が加速されるなら、
    約 70 本の軌道でそれが見える。全数を数えるのは `hot_e` である。
  - 4 ステップに 1 回の記録は、4 倍に圧縮された磁場での電子のジャイロ運動 1 回あたり、約 5、9、14 点
    （L、M、H）になる（本記録の計算）。
  - 計算時間への影響は分かっていない。最初の計算で測る。
- **`spec_up`：** この柱は、B₀ に沿った軸から円筒半径で 204 d_e 以上離れており、電子を 5.8 × 10⁶ 個含む。
  上流の電子温度の時間変化と、数値的な尾の有無が分かる。
- スカラー量（エネルギー、ガウス則、div B）は 37 ステップごとに書く。
- 電子データの解析はまだ書いていない。テーマのフォルダの `analyze.py` として書く。

## 6. 実行前に決めておく確認事項

これは、ある計算を仮説づくりに使えるかどうかを決めるもので、物理の主張ではない。`コード体` の名前は
`examples/blast_wave/analyze.py` の結果の名前である。

1. **エンジンの保存則の契約：** `gauss_max` ≤ 10⁻¹⁰、`divB_max` ≤ 10⁻¹²。
2. **赤道で衝撃波ができている：** 0.3 から 0.7 R_B の間に面のサンプルが 2 個以上あり（`window_samples`）、
   `front_speed` が `fast_speed` を超え、`compression` が 1.5 を超える。これは例の gate の基準である。
   満たさない計算は「衝撃波なし」と報告し、その電子データを衝撃波加速としては読まない。
3. **面が破片より前にある：** `front_debris_distance` > 0。そうでなければ、その計算は、破片が周囲のプラズマを
   突き抜けたものとして扱う。
4. **上流が乱れていない：** `ambient_deviation` ≤ 0.05。
5. **赤道面での対称性：** `asymmetry` ≤ 0.05。
6. **数値加熱：** `spec_up` の電子の平均エネルギーが、ステップ 0 の値の 1.25 倍以内にとどまる。これを超えた
   計算は、β や電子の熱速度に依存する議論には使わない。
7. **粒子配列：** `fill_max` < 0.95。
8. **粒子数：** γ ≥ 1.5 の電子の数を、θ_Bn の 10° の帯ごとに報告する。100 個未満の帯は、そのエネルギーより上の
   スペクトルには使わない。
9. **報告するが基準にはしないもの：**
   - V_d と v_A を単位とした面の速度（その計算の領域）。
   - 破片の最大の広がり。例自身の判定 `passed` にはこの広がりの許容帯が含まれるが、M と H では満たされない
     見込みである。
   - `energy_budget_change`。

## 7. コスト

- **GPU のメモリ：** planner の値は、8 GPU・`precision: compact` で 1 GPU あたり最大 57.9 / 55.9 / 55.3 GB
  （L / M / H）、`fp64` で 63.7 / 61.5 / 60.8 GB である。
  - planner の値は、設定ファイルの診断を数えていない。ヒストグラムは、計算している間、約 7 GB の一時配列を
    必要とする（見積もり）。これは planner が出力用に見込んでいる範囲に収まる。
  - GPU は 1 台 103 GB で（`docs/plans/2026-11-10-bow-shock.md`）、記録にある最大の使用量は 79.8 GB である。
    どちらの精度でも収まる。`compact` は速度のために選んでいる。
- **計算時間：** 8 GPU で 1 本あたり約 15–20 分。これは見積もりで、実測ではない。
  - 例の gate 計算（384³、粒子 9.4 × 10⁸ 個、1 ステップ 73–91 ms）の内訳を、格子 8 倍・粒子 1.9 倍に
    換算した。`compact` で 1 ステップ約 285 ms、`fp64` で約 320 ms になる。
  - そこに出力の分として 10–25 % を足した。
  - 1 セル 2 個しかないので、ここでは粒子より格子のほうがコストが大きい。
- **ディスク：** `hot_e` を除いて、1 本あたり圧縮前で約 13–15 GB。
- **GPU は共用である。** 計算のときは空いていなければならない。オーナーの同意を得て止めるコンテナについては
  `docs/gpu-runbook.md` を見ること。

## 8. あとで実行するときの手順

この段階では行わない。オーナーが計算を始めるときの手順：

1. テーマのブランチの、変更のないコミットの上で作業する。研究ハーネスの実行・記録ツール
   （`research/tools/run.py`、`record_run.py`）は計画だけで、まだ作っていない。それができるまで、計算は
   証拠として記録されない。出力ディレクトリの `run.json` には、コミットと run id が残る。
2. GPU が空いていることを確かめる（`docs/gpu-runbook.md`）。
3. 先に 16 ステップの確認をする。768³ での起動と、ステップ 0 でのすべての診断（空の粒子ダンプを含む）を
   通すためである。そのダンプが失敗したら `--set output.codec=none` を付ける。
4. 最初に L を走らせる。例の gate の領域に最も近く、計算時間、メモリ、数値加熱、`hot_e` のサイズの最初の
   実測が得られる。次に H、最後に M。

```bash
# 16 ステップの確認
python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml \
  --set overrides.n_steps=16 --out output/research/spherical-shock-injection/check-L
# 計算と、例の解析
python -m pypic.launch -n 8 examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml \
  --out output/research/spherical-shock-injection/explore-L
python examples/blast_wave/analyze.py output/research/spherical-shock-injection/explore-L
```

計画の表示だけなら、GPU のない計算機でもできる（何も書き出さない。numpy、pyyaml、scipy が必要）：

```bash
PYPIC_WORLD_SIZE=8 PYTHONPATH=. python examples/blast_wave/run.py \
  --config research/topics/spherical-shock-injection/configs/explore-L.yaml --plan-only
```

planner が確かめるのは、診断の名前とキーの名前である。ヒストグラムの設定と `select` のキーは、計算の開始時に
初めて確かめられる。ここでは `pypic/output/diagnostics.py` と手で突き合わせて確かめた。

## 9. 今回外したものと、それに必要なもの

- **速度を保つ面（ピストン駆動）。** 破片の質量をずっと大きくすれば、面のマッハ数を保てる。planner は箱を
  R_B に結び付けているので、箱が大きくなりすぎる。テーマ独自の planner が必要になる。これは 4.3 節の最初の
  リスクへの直接の答えである。
- **M_A が 20–30。** 磁場を面内に置いた、薄い準 2 次元の箱なら可能である：円筒形の面で、すべての θ_Bn を
  持つ。これも独自の planner が必要で、球ではなくなる。
- **冷たい電子**（残骸の β と、電子の熱速度に対するマッハ数）：0.2 d_e 以下の格子が必要で、この大きさの
  3 次元の箱では無理である。
- **別の質量比と衝撃波速度。** スケーリングについて何かを言う前に必要である（調査 7 節）。
- **同じ局所パラメータでの平面の計算。** 曲がった面が何を加えるかを切り分けるため。
- **必要な場所に置くテスト電子**（`tracks_e` の間引きの代わり）、半径と角度でビンに分けたスペクトル、領域と
  エネルギーを組み合わせた選別。

## 10. 3 本の設定に対する planner の出力

2026-11-12 に、8 節の最後のコマンドで出力した（CPU のみ。計画は 8 GPU 用）。

**L**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.285714 along z, beta=0.49, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0714 c, fast speed 0.0848 c
debris (on): M_A=3.5 (fast Mach 2.95), V_d=0.2500 c, R_m=64 d_e = 16 d_i, R_B=147.5 d_e, gyroradius 14 (0.219 R_m); shell at 21.12, width 6.40, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 67.2 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.25 R_B/V_d = 13.2 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 74,672,289 of the debris (block of 8,998,912 cells), capacity factor 1.5; memory up to 57.9 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 13 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11   23.3 MB    257 MB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~356 GB
```

**M**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.166667 along z, beta=1.44, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0417 c, fast speed 0.0618 c
debris (on): M_A=6 (fast Mach 4.05), V_d=0.2500 c, R_m=44.8 d_e = 11.2 d_i, R_B=147.9 d_e, gyroradius 24 (0.536 R_m); shell at 14.78, width 4.48, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 51.9 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.25 R_B/V_d = 7.7 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 25,607,222 of the debris (block of 3,112,136 cells), capacity factor 1.5; memory up to 55.9 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 9 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11   68.7 MB    756 MB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~356 GB
```

**H**

```text
=== pypic blast-wave planner (debris shell around a magnetised body) ===
ambient: m_i/m_e=16, B0=0.114943 along z, beta=3.03, v_te=0.1 c (Debye length 0.20 dx), v_A=0.0287 c, fast speed 0.0539 c
debris (on): M_A=8.7 (fast Mach 4.64), V_d=0.2500 c, R_m=34.8 d_e = 8.7 d_i, R_B=147.2 d_e, gyroradius 34.8 (1.000 R_m); shell at 11.48, width 3.48, peak density 16.5 in 4 beam(s) per species (omega_pe dt <= 1.05, 0.00 % inside the body); extent after Le et al. 2021, Eq. 4 with D = 3: 43.5 d_e
body: OFF (no sphere, no dipole)
grid (768, 768, 768), dx=0.5, L=384 d_e (half-width 1.30 R_B), dt=0.25, 2960 steps = 1.26 R_B/V_d = 5.3 / Omega_ci
1,811,939,328 ambient particles (2 per cell and species), about 12,002,231 of the debris (block of 1,481,544 cells), capacity factor 1.5; memory up to 55.3 GB/GPU on 8 GPU(s), rank grid (1, 2, 4); ids 6.0e-03 of 2^53
faces: six open faces with flux injection of the ambient species; filter 4+comp, Marder zone 0.125, npml 12
output: 3-D moments in blocks of 7 cells every 296 steps, slices z and y through the centre every 74 steps, history every 37 steps
output plan (uncompressed sizes; zstd typically saves 10-50 % more):
  diagnostic       kind         every    dt_out per period outputs      each     total  rule
  coarse           fields         296        74        0.1      11    144 MB   1.58 GB  296 steps
  slice_z          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  slice_y          slice           74      18.5        0.3      41   47.2 MB   1.93 GB  74 steps
  spec_meridional  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_equatorial  histogram      148        37        0.2      21   56.6 MB   1.19 GB  148 steps
  spec_all         histogram       74      18.5        0.3      41   3.58 kB    147 kB  74 steps
  spec_up          histogram       74      18.5        0.3      41   1.79 kB   73.5 kB  74 steps
  hot_e            particles      148        37        0.2      21   16.3 GB    342 GB  148 steps
  tracks_e         tracks           4         1        6.3     741   9.19 MB   6.81 GB  4 steps
  scalars every 37 steps; total ~357 GB
```
