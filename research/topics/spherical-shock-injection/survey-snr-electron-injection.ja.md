# 超新星残骸の衝撃波における電子注入：先行研究の調査

テーマ `spherical-shock-injection`。2026-11-12 作成。状態：最初の調査で、レビュー前。

1 ページの要約は `survey-snr-electron-injection-summary.ja.md` にある。
これは `survey-snr-electron-injection.md`（英語、原本）の日本語版である。節番号と内容は英語版に対応する。
食い違いがあれば英語版に従う。

**対象。** 超新星残骸（SNR）の衝撃波における電子注入問題を扱う。観測が何を要求するか、どんな機構が提案されて
いるか、運動論的シミュレーションが何を示したか、何が未解決かをまとめる。太陽圏の衝撃波は、SNR の衝撃波の
検証材料になる場合だけ扱う。

**読み方。**
- `park2015simultaneous` のような名前は citekey である。ノートは `research/literature/<citekey>.md` にあり、
  書誌情報は 10.1 節にある。論文中の位置（Sect.、Fig.、Eq.、Table）は、読んだ arXiv 版のものを英語表記のまま示す。
- 「未読」は、この調査では読んでいない文献を指す。読んだ文献を通して引用し、書誌情報はその文献の記載どおりに
  10.2 節に載せる。
- 「本調査の計算」は、論文の数値からここで導いた値である。論文の記述ではない。
- 略語：θ_Bn は衝撃波法線と上流磁場のなす角。M_A と M_s はアルフベン・マッハ数と音速マッハ数。SSA は
  サーフィン加速、SDA は衝撃波ドリフト加速、SSDA は確率的衝撃波ドリフト加速、DSA は拡散衝撃波加速。
  K_ep は加速粒子の電子・陽子比。PIC は particle-in-cell。

## 1. 要約

1. **問題。** DSA が加速できるのは、衝撃波を横切れる粒子だけである。実際の SNR 衝撃波では、熱的電子とその
   閾値の間にエネルギーで 10⁵–10⁶ 倍の隔たりがある（`bohdan2023electron` Eq. 1）。電子がこの隔たりをどう
   越えるかが注入問題である。閾値が何であるかについて、文献は一致していない（2 節）。
2. **観測**が与えるのは注入の結果であって、過程ではない。若い SNR には 10–100 TeV の電子がある
   （`vink2012supernova` Sect. 6.3.1）。Tycho と SN 1006 のモデルは K_ep ≈ 10⁻³ を必要とする
   （`morlino2012strong` Sect. 3、`mcginness2026sn` Table 1）。注入された電子の割合を測った観測はない。
3. **残骸での角度依存は決着していない。** SN 1006 については、多くの研究が効率のよい加速を準平行の領域
   （極冠）に置く。一方、双極型の電波残骸の集団は準垂直での加速を支持する（`west2017connection` Sect. 5）。
   残骸全体のモデルはどれも、効率の角度依存 ε(θ_Bn) を仮定として入れている（3.3 節）。
4. **機構は角度によって違う**（4 節）：
   - 準垂直：Buneman 波でのサーフィン、散乱を伴う衝撃波ドリフト（SSDA）、ホイッスラー波、Weibel 乱流中の
     リコネクションと確率的散乱。
   - 斜め：電子フォアショックと、そこでの電子音波・ホイッスラー波。
   - 準平行：ミラー反射のあと、衝撃波と陽子が駆動する波の間での散乱。

   DSA と整合する電子スペクトルに達したのは、1 次元（`park2015simultaneous`、`xu2020electron`、
   `kumar2021nonthermal`、`gupta2024electron`）、薄い 2 次元の箱（`arbutina2021nonlinear`）、弱相対論的な速度での
   2 次元（`crumley2019kinetic`）である。このうち、2 次元以上での準垂直衝撃波は一つもない。
5. **一組のパラメータで θ_Bn の全域を扱った運動論的シミュレーションはない。** また、SNR 相当のパラメータで、
   大きなスケールで曲がった衝撃波面を扱ったものも見つからなかった（6 節）。曲がった面は地球のバウショックに
   ついてだけ扱われており、その文献はここでは読めなかった。
6. **数値上の制約が答えを左右する**（7 節）：小さくした質量比、0.03–0.4 c の衝撃波速度、1 次元や 2 次元の
   幾何。2 次元では、磁場の向きによって、どの不安定性が存在できるかが決まる。機構ごとの閾値は、質量比の
   異なるべきでスケールする。
7. **未解決**（8.2 節、8.3 節）：ε(θ_Bn) の形。現実の質量比と衝撃波速度でどの機構が支配するか。準垂直と
   準平行のどちらの領域が電子を多く注入するか。2 次元以上で、前段加速から DSA へどう移るか。

## 2. 注入問題を数値で見る

**「注入された」の意味。** 文献によって閾値が違う：

| 閾値 | 内容 | 出典 |
|---|---|---|
| ジャイロ半径 ≥ 衝撃波の厚み | 実際の衝撃波ではエネルギーで 10⁵–10⁶ 倍の隔たり。斜め衝撃波の PIC 計算では 10³ 倍 | `bohdan2023electron` Eq. 1, Sect. 4.3 |
| 運動量 | p_e ≫ m_i V_sh。V_sh = 0.01 c で約 10 MeV/c | `pohl2020pic` Eq. 34、`bohdan2019kineticb` Sect. 4 |
| 拡散長 ≲ 衝撃波の厚み | ジャイロ半径の条件は「最も悲観的」。注入は 0.1–1 MeV | `amano2022nonthermal` Sect. 2.5, Eq. 3 |
| 陽子の注入運動量 | 電子は陽子の p_inj ≈ 3 m_p v_u の近くで DSA に入る | `park2015simultaneous` Fig. 3 |
| 速さ | 電子は、速さが v_sh を超えれば、ジャイロ半径によらず注入される | `gupta2025speed` Sect. II B, Fig. 2 |
| 上流構造のスケール | 高 M_A では、ジャイロ半径 ≈ 急峻化した波（SLAMS）の λ/2π | `zekovic2025slams` Eq. 1 |
| 熱的運動量の倍数 | p_inj = ξ p_th、ξ ≈ 3.2–3.7 で両粒子種に共通 | `arbutina2021nonlinear` Eqs. 4, 12, Sect. IV |

**SNR 衝撃波のパラメータ。** シミュレーションが合わせるべき値：
- 衝撃波速度は、X 線シンクロトロンの縁を持つ残骸で 2000–6000 km/s（`vink2012supernova` Sect. 11.3）。
  これは 0.007–0.02 c にあたる（本調査の計算）。SN 1006 は 5000 ± 400 km/s、北西では 3000 ± 400 km/s
  （`mcginness2026sn` Sect. 2.1）。Tycho はモデルで 4990 km/s（`morlino2012strong` Sect. 2.1）。
- 磁場は、星間空間で 3–5 μG（`reynoso2013radio` Sect. 4、`vink2012supernova` Sect. 11.3）、縁で 100–600 μG
  （`vink2012supernova` Sect. 11.2）。
- 密度は 0.02–0.3 cm⁻³（`mcginness2026sn` Sect. 2.1、`morlino2012strong` Sect. 2.1）。
- マッハ数を与えているのは `morlino2012strong` だけで、増幅後の磁場で M_s ≃ 443、M_A ≃ 15（Sect. 2.3）。
  増幅前の 5 μG と 0.3 cm⁻³ を使うと、遠方上流で M_A ≈ 250、ω_pe/Ω_ce ≈ 350、β ≈ 0.8 になる（本調査の計算）。
- 衝撃波が超光速になり、電子が磁力線に沿って衝撃波を追い越せなくなるのは、θ_Bn > 89° だけである
  （`amano2022nonthermal` Sect. 2.1）。
- 半径と衝撃波の厚みの比は約 10¹²（`vink2012supernova` Sect. 5.2 からの本調査の計算）。

**シミュレーションのパラメータ**（5 節）：質量比 m_i/m_e は多くが 16–400。衝撃波速度は 0.03–0.4 c。M_A は
多くが 5–70 で、2 次元では 150 まで（`bohdan2021magnetic`）、1 次元では 300 まで（`zekovic2025slams`）。
0.26 c では、θ_Bn ≈ 75° より上で衝撃波が超光速になる（`morris2022preacceleration` Sect. 2）。一様磁場中の
球面のうち超光速になる部分の面積比は v_sh/c で、0.26 c なら 26 %、5000 km/s なら 2 % 未満である
（本調査の計算）。

## 3. 観測からの制約

### 3.1 電子は SNR の衝撃波で加速されている
- 薄い X 線シンクロトロンの縁は、若い残骸の前進衝撃波に 10–100 TeV の電子があることを示す
  （`vink2012supernova` Sect. 6.3.1）。最初のシンクロトロン殻は SN 1006 で見つかった（Koyama et al. 1995、未読。
  `vink2012supernova` と `reynolds2004microphysics` が引用）。
- カットオフの光子エネルギーは hν = 1.4 η⁻¹ (V_s / 5000 km/s)² keV で、η は Bohm 値を単位とする拡散係数
  である（`vink2012supernova` Eq. 52）。
- 熱的電子は冷たいままである：若い残骸すべてで kT_e ≤ 5 keV。SN 1006 では kT_e ≈ 1.5 keV に対し、O VII は
  528 ± 150 keV（`vink2012supernova` Sect. 11.1）。
- `vink2012supernova`（Sect. 5.5）は注入を "outstanding uncertainty" と呼んでいる。

### 3.2 電子はどれだけあるか
- Tycho：電波と X 線のシンクロトロン放射から K_ep = 1.6 × 10⁻³。これはモデルの自由な規格化である
  （`morlino2012strong` Sect. 3）。同論文は、電子注入の理論は "still lacking" と述べている。
- SN 1006：象限ごとのモデルで、4 象限すべて K_ep = 9.3 × 10⁻⁴（`mcginness2026sn` Table 1）。このモデルでは
  電子は作り方の上で陽子に結び付けられている。同論文は、文献中の K_ep の範囲を 10⁻²–10⁻⁴ としている
  （Sect. 2.2）。
- 地球での値は K_ep ≈ 0.01 で、若い SNR は 10⁻³ 以下を示唆する（`park2015simultaneous`、Introduction）。
- これらのフィットでは、K_ep は磁場と縮退している（`morlino2012strong` についての本調査の解釈）。
- 熱的電子のうち注入される割合を測った論文はない。

### 3.3 殻のどこで起きるか：角度依存
一様磁場中の球形の残骸には、すべての θ_Bn がある。SN 1006 には明るい縁が二つあり、それが極冠（準平行）なのか
赤道帯（準垂直）なのかが問題になる。

**モデルへの角度依存の入れ方。** 次の法則はどれも運動論から導かれたものではない：

| 法則 | 形 | 出典 |
|---|---|---|
| 準平行 / 準垂直 / 等方 | ε ∝ cos²Θ、sin²Θ、定数。Θ は衝撃波後面の磁場に対してとる | `bocchino2011constraints` Sect. 2、`west2017connection` Sect. 2。起源は Fulbright & Reynolds 1990、Leckband et al. 1989（未読） |
| 加速率に角度依存を入れる | κ = κ_∥ [cos²θ_Bn + sin²θ_Bn / (1 + η²)]。垂直衝撃波のほうが速い | `reynolds2004microphysics` Eq. 3、Jokipii 1987（未読） |
| 陽子の階段関数 | ζ_p = (0.15/2) [tanh((θ_crit − θ)/δ) + 1]、θ_crit = 45°、δ = 10° | `winner2020evolution` Eq. 1。`caprioli2014simulations` に基づく |
| 電子の階段関数、3 通りの試行 | 同じ tanh で ζ_∥ と ζ_⊥ をつなぐ。ζ_∥/ζ_⊥ = 10、1、0.1 | `winner2020evolution` Eq. 4, Fig. 1 |
| 加速粒子の圧力 | tanh の階段。45°（熱的陽子、12 %）と 70°（再加速された種粒子、6 %） | `giuffrida2022supernova` Methods, Eq. 1 |
| 象限ごとのフィット | 効率は極冠で 20 %、それ以外で 0.4–0.8 %。K_ep は一定 | `mcginness2026sn` Table 1 |

**証拠。**

| 文献 | データまたは方法 | 結果 | その文献自身の留保 |
|---|---|---|---|
| Rothenflug et al. 2004（未読） | SN 1006 の X 線の中心対縁の比 | 赤道帯を排除する | `bocchino2011constraints` Sect. 1 経由 |
| `reynolds2004microphysics` | 角度依存の拡散を入れた合成画像 | 明るい縁は "almost certainly perpendicular"（Conclusions） | 平均磁場 10 μG を仮定 |
| `bocchino2011constraints` | 3 次元 MHD と SN 1006 の電波画像 | 磁場が視線に対し 38° ± 4° で、磁場強度に勾配があるとき、準平行が最もよく合う（Sect. 3.4.1） | 一様磁場では準平行は合わない（Sect. 3.3）。準垂直は "worse" だが排除されない（Sect. 3.4.2） |
| `reynoso2013radio` | SN 1006 の電波偏波 | 明るいローブ：磁場は動径方向で乱れており、偏波率 0.17。南東の縁：接線方向の整った磁場で 0.6 ± 0.2（Sects. 3.1, 3.2） | 勾配モデルを "somewhat contrived" と呼ぶ（Sect. 1） |
| `zhou2023magnetic` | 北東の縁の X 線偏波 | 下流の磁場は動径方向。偏波度 22.4 ± 3.5 %（Table 1） | 準平行と "consistent with"。動径磁場と一様磁場は区別できない（Sect. 3.2） |
| `giuffrida2022supernova` | 衝撃波後面の密度の方位角依存 | 圧縮比は θ ≈ 0° に向かって 4 から約 7 に上がる（Fig. 3b） | ハドロンについての結果。電子の数はない |
| `winner2020evolution` | 電子スペクトルを入れた 3 次元 MHD | 画像を再現するには、電子の効率が準平行で 10 倍以上高い必要がある（Sects. 3.1, 7） | 試した法則は 3 通り。磁場モデルは外から与えている |
| `west2017connection` | 銀河磁場モデル中の双極型電波残骸 32 個 | 準垂直の重みが約 75 % と合う（Sect. 5）。準平行は少なくとも 10 個で合わない（Sect. 3） | 目視による比較。SN 1006 は準平行のほうが合う（Sect. 4.2） |

**食い違い。**
- SN 1006 は、多くの研究が準平行とみなす。`reynolds2004microphysics` とそれが従う研究は準垂直とみなす。
- 双極型残骸の集団は逆を指す（`west2017connection`）。双極型残骸の多くは軸が銀河面に沿っており、SN 1006 は
  "rather peculiar" である（`bocchino2011constraints` 脚注 5）。
- `winner2020evolution`（Sect. 7）は、自らの結果が、電子について準垂直衝撃波を支持する 1 次元 PIC 計算の
  外挿を許さない（"precludes extrapolation"）と述べる。同論文の準垂直の試行は、まさにその計算（`xu2020electron`）を
  動機としており、準平行の側は `park2015simultaneous` である（Sect. 2.2）。
- SN 1006 の磁場はモデルによって違う：準垂直の縁で 140 μG（`winner2020evolution`）に対し 6 μG
  （`mcginness2026sn`）。

**論文自身が未解明と述べていること。** 電子の法則の形は第一原理からは分かっていない
（`winner2020evolution` Sect. 2.2）。電子注入は "is not completely understood" であり、準垂直領域での機構は
未解決である（`mcginness2026sn` Sects. 2.2, 3.2）。

## 4. 機構と、運動論的シミュレーションが示すこと

`bohdan2023electron`（Sects. 2–5）は、分野を角度で整理している：垂直（75–90°）、斜め（50–75°）、
平行（50° 未満）。この節はその順に従う。

### 4.1 高マッハ数の準垂直衝撃波

**Buneman 波でのサーフィン（SSA）。** 反射されたイオンがフット領域で流入電子に対してドリフトし、Buneman
不安定性を駆動する。波に捕捉された電子が運動電場からエネルギーを得る。
- 起源：Shimada & Hoshino 2000、McClements et al. 2001（未読）、`hoshino2002nonthermal`。
- 条件：
  - 不安定性：M_A ≥ [(1+α)/2] (β_e m_i/m_e)^(1/2)。α ≈ 0.2 は反射イオンの割合
    （`matsumoto2012electron` Eq. 5）。
  - 捕捉：M_A ≥ (1+α) (m_i/m_e)^(2/3)（`matsumoto2012electron` Eq. 8）。m_i/m_e = 100 で約 26、1836 で
    約 180 になる（`wieland2016nonrelativistic` Sect. 2.2）。
  - それ以前の 1 次元の条件は、現実の質量比で M_A > 43–75 だった（`hoshino2002nonthermal` Eq. 8）。
- 次元によって結果が変わる：
  - 1 次元：平面波が電子を捕捉する（`amano2007electron`）。
  - 磁場が面外の 2 次元：サーフィンは確率的な過程として働く（`amano2009electron`、
    `matsumoto2012electron`、`matsumoto2013electron`）。
  - 磁場が面内の 2 次元：Buneman 波は抑えられる。前段加速される割合は約 0.45 % で、面外では 6.8 %
    （`bohdan2019kinetic` Table 2）。
  - 3 次元：サーフィンは存在し、第二段階に粒子を渡す（`matsumoto2017electron`）。
- 否定的な結果：
  - 現実の質量比での 2 次元静電計算では加熱しか起きない。波のポテンシャルが 1 波長ぶんしか揃わないため
    である（`ohira2007absence` Sects. 2.4, 3）。
  - 最大エネルギーは衝撃波速度で決まり、質量比によらない。そのため 1836 では無視できる
    （`riquelme2011electron` Appendix）。
  - 電子 1 個あたりの平均の利得は質量比によらないが、下流の温度は質量比とともに上がる。そのため、非熱的な
    尾のうちサーフィンによる部分は、m_i/m_e が 50 から 400 になると、電子の 0.15 % から 0.0015 % に下がる
    （`bohdan2019kineticb` Sect. 3.3）。

**衝撃波ドリフト（SDA）と確率的衝撃波ドリフト（SSDA）。** 衝撃波でミラー反射された電子は、運動電場から
エネルギーを得る。衝撃波層の中でピッチ角散乱があると、この過程が繰り返される。
- 起源：Wu 1984、Leroy & Mangeney 1984（未読）。
- 1 回の反射で Δv = 2 u₀ / cos θ_Bn（`katou2019theory` Eq. 1）。θ_Bn が大きいほど利得は増え、反射される
  電子の数は減る（`amano2007electron` Sect. 1）。
- 散乱がある場合：べき分布 N(v) ∝ v^(−p)、p = 1 + 3L/L_sh ≈ 4（`katou2019theory` Eq. 20）。最大エネルギーは
  ε_max ~ ½ m_i u_sh² D_μμ/Ω_ce で、u_sh = u₀/cos θ_Bn（Eq. 21）。地球のバウショックで約 11 keV（Eq. 34）、
  その 10 倍の速度で約 1 MeV（Sect. 4）になる。
- その場観測による裏付け：M_A ≃ 8.9、θ_Bn ≃ 85° のバウショック通過 1 例で、カットオフが 22 ± 1 keV。予測は
  約 20 keV（`amano2020observational` p. 4, Fig. 4）。
- 閾値：
  - 反射電子が駆動するホイッスラー波：M_A ≳ (cos θ_Bn / 2) (β_e m_i/m_e)^(1/2)（`amano2010critical` Eq. 2）。
  - SSDA を通じた注入：M_A / cos θ_Bn ≳ 55。そこで仮定した散乱率と厚みでの値（`amano2022theory` Eq. 14）。
  - 反射電子による波の励起：M_A / cos θ_Bn ≳ (α/4) (m_i/m_e)^(1/2) ≈ 50（`wang2026generation` Eq. 22）。
  - M_A / cos θ_Bn が非常に大きいと、反射は再び弱まる：反射電子の密度は
    exp[−(M_A/cos θ_Bn)² (m_e/m_i) / ((B₂/B₁ − 1) β_e)] で下がる（`wang2026generation` Eq. 24）。
- θ_Bn = 74.3° の 3 次元 PIC：まずサーフィン、次に Weibel 乱流による散乱を伴う SSDA。指数 3.5 のべき分布。
  γ > 5 の電子は上流へ逃げる（`matsumoto2017electron` Figs. 3, 4）。
- 質量比を小さくすると、この注入は人工的に起きやすくなる（`amano2022theory` Eq. 21, Sect. 4.4）。

**ホイッスラー波。** フット領域の斜めホイッスラー波が電子を磁力線方向に加速する。電子は上流へ逃げ、
また戻ってくる（`riquelme2011electron`）。
- 条件：M_A ≲ (m_i/m_e)^(1/2)。m_i/m_e = 1600 まで調べられている（Sect. 3.1, Fig. 12）。
- 結果：電子の数で 2–10 %、エネルギーで 10–30 % が、指数 2.7–4 の尾に入る（Sect. 4.4）。
- SNR では M_A ≲ 20、すなわち上流磁場が増幅されていることが必要になる（Sect. 5）。
- これは、M_A が (m_i/m_e)^(2/3) より大きいことを求めるサーフィンの条件と矛盾する。

**Weibel 乱流：リコネクションと確率的散乱。** M_A が数十以上では、イオン Weibel 不安定性が衝撃波の構造を
決める。
- M_A ~ 130 で磁場が面内の垂直衝撃波は Weibel 媒介で、電子加速は見られない（`kato2010nonrelativistic`）。
- Weibel 電流シートのリコネクション（Matsumoto et al. 2015、未読）が、m_i/m_e ≥ 100 では尾の大部分を担う：
  非熱的な割合 0.55–0.7 % のうち電子の 0.4–0.55 %（`bohdan2020kinetic` Table 3）。残りは確率的な Fermi 型の
  散乱によるとされる（`bohdan2019kineticb`）。
- 垂直衝撃波での非熱的な割合は約 0.5 % で、到達するのは高々 10–100 kT である（`bohdan2023electron`
  Sect. 3.3）。これは前段加速であり、注入の閾値には遠い。
- Weibel 磁場は B_sh ≈ 2 √M_A B₀ で成長する。これにより注入運動量が √M_A で上がる
  （`bohdan2021magnetic` Eqs. 2, 6）。

**1 次元での電子駆動の波と DSA。**
- θ_Bn = 63°、M_A = 63、m_i/m_e = 100：反射電子が遠方上流まで波を駆動し、p⁻⁴ のスペクトルで DSA に達する。
  電子の数で約 7 %、エネルギーで約 20 %。陽子は注入されない（`xu2020electron` Sect. 3, Fig. 4）。
- θ_Bn = 75°、M_A = 63、m_i/m_e = 256：両粒子種が DSA を始める。反射は衝撃波の再形成の位相に依存する
  （`kumar2021nonthermal` Sects. 3, 4）。
- どちらも 1 次元である。3 Ω_ci⁻¹ の 3 次元での確認では、反射される電子が "by a factor of few" 少ない
  （`xu2020electron` Fig. 6）。

### 4.2 斜め衝撃波：電子フォアショック
- 反射には、衝撃波静止系で v_e cos θ_Bn ≥ v_sh が必要である（`morris2022preacceleration` Sect. 4.1.1）。
- 反射される割合は角度とともに下がる：30° で約 5 %、63° で 0.05 %、74.3° でゼロ（m_i/m_e = 64、
  M_A = 21、v_sh = 0.264 c。`morris2022preacceleration` Sect. 4.1, Figs. 4, 5）。63° で磁場を面外にすると 4 %。
- 外側のフォアショック：電子音波。θ_Bn ≤ 45° で、上流電子の約 1 % を衝撃波に届かなくする（Sect. 6）。
  効果は小さい。
- 内側のフォアショック：斜めホイッスラー波が成長して、イオン慣性長の 1–5 倍の空洞になる。空洞は反射電子を
  閉じ込めるが、エネルギーの利得は衝撃波で起きる。下流では、電子の 0.8 % がエネルギーの 7 % を持ち、
  指数は −2.8 まで（`morris2023preacceleration` Sects. 4.6, 4.8）。DSA には達していない。
- フォアショックが定常になるには 125–300 Ω_ci⁻¹ かかる（`morris2023preacceleration` Sect. 4.5、
  `morris2022preacceleration` Sect. 5）。5 節の 2 次元計算の長さは 6–52 Ω_ci⁻¹ である。
- 上流にあらかじめ乱流がある場合：垂直衝撃波では変化がない（`fulat2023kinetic` Table 1, Fig. 9）。60° では、
  非熱的な割合が数で 0.14 % から 0.23 % に上がり、傾きは同じである（`fulat2026pic` Sect. 3.4）。これは
  1 組の計算による。

### 4.3 準平行衝撃波
- **基準となる描像**（`park2015simultaneous`）：1 次元、θ_Bn = 30°、M_A = 20、m_i/m_e = 100。
  - 電子はミラー反射され、衝撃波ドリフトでエネルギーを得る。
  - その後、衝撃波と、陽子が駆動する Bell 波との間に捕捉され、陽子の注入運動量に達する。
  - そこから両粒子種が p⁻⁴ に従う。
  - K_ep = 3.8 × 10⁻³。m_i/m_e = 400 では 5.5 × 10⁻³。
- **同じ角度での別の結果**（`kato2015particle`）：電子の注入は 1 回のバーストだけで、エネルギーの指数は
  約 3.0（Sect. 3.3, Figs. 15, 16）。
- **30° での 1 次元の系統的計算**（`gupta2024electron`）：
  - 非熱的電子の割合は数で高々 0.5 %。
  - スペクトルが p⁻⁴ に達するのは M_A, M_s ≳ 10 のときだけ（Sects. 5, 6）。
  - この割合は質量比によらない（Fig. 7d）。
  - 現実の質量比での電子・イオン比は高々 10⁻²。これは式による値である（Eq. 7）。
- **速さの閾値**（`gupta2025speed`）：電子は、速さが v_sh を超えると DSA に入る。反射には、すぐ上流での
  電子の実効的な音速マッハ数が約 3 未満であることが必要（Eq. 4）。約 10 % が反射され、そのうち注入されるのは
  高々半分（App. C）。K_ep ~ 10⁻³（Eq. 5）。
- **平行衝撃波での競合する機構**（`shalaby2022mechanism`）：中間スケール不安定性のイオンサイクロトロン波が
  電子を散乱する。衝撃波面での M_A が (m_i/m_e)^(1/2)/2 未満であることが必要（Eq. 2）。そのため質量比を
  小さくするとこの機構は消える：m_i/m_e = 100 では 2–3 桁小さくなる（Sect. 6）。
- **高 M_A**（`zekovic2025slams`）：陽子が駆動する波が急峻化して大振幅の構造（SLAMS）になる。B_⊥/B₀ は
  M_A = 20 で 3–5、M_A = 300 で 30 以上（Table 1）。これが電子を捕捉し、より速い Fermi 型の過程が、DSA の
  閾値まで、傾き 5–5.7 の急なスペクトルを作る（Eq. 7, Fig. 4）。DSA への移行に達した計算はない（Sect. 5）。
- **2 次元での領域の切り替わり**（`jikei2026dependence`）：θ_Bn = 20° で、M_A が約 100 未満の衝撃波は
  Bell 型の磁場が支配し、衝撃波エネルギーの約 1.5 % が非熱的電子に入る。それより上では Weibel 支配になり、
  約 15 % になる（Sect. 3.3）。計算は準相対論的な速度で行われている。
- **モデル用の処方**（`arbutina2021nonlinear`）：電子と陽子に同じ ξ を使うと、5000 km/s で
  K_ep = 0.0007–0.0019 になる（Sect. IV）。
- **弱相対論的、2 次元**（`crumley2019kinetic`）：10° では p^(−4.2) の電子の尾ができ、エネルギーの
  5 × 10⁻⁴ を持つ。55° はその速度では超光速で、尾はできない（Sects. 4, 6）。

### 4.4 機構の連鎖
高エネルギーに達した研究は、どれも連鎖になっている：電子を熱的な集団から引き上げる第一段階と、それを
引き継ぐ第二段階である。

| 文献 | 第一段階 | 第二段階 | 到達点 |
|---|---|---|---|
| `amano2007electron` | サーフィン | ミラー反射と衝撃波ドリフト | γ ~ 2、上流へのビーム |
| `matsumoto2012electron` | サーフィン（複数回） | オーバーシュートでの断熱的な利得 | γ ~ 7–9 |
| `matsumoto2017electron` | サーフィン | Weibel 散乱を伴う SSDA | γ − 1 ~ 5、上流への脱出 |
| `bohdan2020kinetic` | （サーフィンは小さい） | リコネクション、確率的散乱 | 非熱的 0.5 % |
| `riquelme2011electron` | ホイッスラー波の平行電場 | 脱出、帰還、ドリフト型の利得 | Γ ≈ 50 |
| `morris2023preacceleration` | 反射 | ホイッスラー空洞での閉じ込め、衝撃波での利得 | γ ≈ 13 |
| `xu2020electron` | フットでの加熱、反射 | 電子駆動の波、衝撃波ドリフトの繰り返し | DSA |
| `park2015simultaneous` | 反射、衝撃波ドリフト | 陽子駆動の波による捕捉 | DSA |
| `zekovic2025slams` | 反射 | 衝撃波と SLAMS の間の Fermi 過程 | DSA の閾値の手前 |

### 4.5 衝撃波の微視的物理によらない注入
種となる電子を原子過程や核過程から供給する提案が二つある：加速された原子核の電離（Morlino 2009）と、
イジェクタ中の放射性崩壊（Zirakashvili & Aharonian 2011）。検索で見つけたが読んでいない（10.3 節）。
このテーマの範囲外である。

## 5. 運動論的シミュレーションの一覧

断りのないものは PIC である。「面内」「面外」は、2 次元での磁場の向きを示す。マッハ数は論文の記載どおりで、
基準系は論文によって違う（ノートを参照）。結果は電子についてのものである。

**垂直・準垂直（θ_Bn ≥ 60°）**

| 文献 | 次元 | m_i/m_e | M_A | β | θ_Bn | v_sh/c | 結果 |
|---|---|---|---|---|---|---|---|
| `amano2007electron` | 1 次元 | 100 | ≈15 | 各 0.08 | 70–85° | ≈0.075 | 反射率は 70、75、80° で ≈5e-2、1.5e-2、2e-3。85° ではゼロ（Fig. 7, Sect. 3.3） |
| `amano2009electron` | 2 次元、面外 | 25 | ≈14 | 各 0.5 | 90° | 0.2（箱の系） | 指数 2.0–2.5 の尾（Fig. 4） |
| `ohira2007absence` | 2 次元静電、周期的なフットのモデル | 1836 | 記載なし | 記載なし | 90° | 0.02 | 加熱のみ |
| `matsumoto2012electron` | 2 次元、面外 | 25, 100 | 14–30 | β_e 0.5, 4.5 | 90° | 0.2（箱の系） | 捕捉条件と低い β_e を満たすと尾ができる（Table 1） |
| `matsumoto2013electron` | 2 次元、面外 | 225 | 44.8 | β_e 0.5 | 90° | 記載なし | γ_max ~ 12（Fig. 2b） |
| `matsumoto2017electron` | 3 次元 | 64 | 20.8 | 各 0.5 | 74.3° | 0.26 | 指数 3.5 のべき分布（Fig. 4d） |
| `kato2010nonrelativistic` | 2 次元、面内 | 30 | ~130 | ~26 | 90° | 0.33 | マクスウェル分布。T_e/T_i ≈ 0.33（Fig. 14） |
| `riquelme2011electron` | 2 次元の面内と面外、薄い 3 次元 | 25–1600 | 3.5–14 | 各 0.005–0.5 | 45–90° | 0.042, 0.14 | 数で 2–10 %。指数 2.7–4。90° では何も出ない（Sect. 4.4） |
| `wieland2016nonrelativistic` | 2 次元、磁場は面に対し 45° | 50 | 28 | β_e ≪ 1 | 90° | 0.13, 0.39 | 尾はない（Figs. 16, 17） |
| `bohdan2017electron` | 2 次元、面内・45°・面外 | 100 | 32–36 | 5e-4, 0.5 | 90° | 0.26–0.29 | 非熱的 0.2–0.7 %（面内、45°）、4–7 %（面外）（Table 2） |
| `bohdan2019kinetic`、`bohdan2019kineticb`、`bohdan2020kinetic` | 2 次元、面内（1 本は面外） | 50–400 | 23–69 | 5e-4, 0.5 | 90° | 0.26 | 非熱的 0.1–0.7 %。m_i/m_e ≥ 100 では大部分がリコネクションから |
| `bohdan2021magnetic` | 2 次元、面内 | 50–400 | 23–150 | 5e-4, 0.5 | 90° | 0.26 | 磁場増幅のみ |
| `fulat2023kinetic` | 2 次元、面内、上流に乱流 | 100 | 32 | 1 | 90° | 0.264 | 乱流による変化なし |
| `fulat2026pic` | 2 次元、面外、上流に乱流 | 100 | 32 | 1 | 60° | 0.264 | 乱流で 0.14 % → 0.23 %（Sect. 3.4） |
| `morris2023preacceleration` | 2 次元、面外 | 50 | 30 | 1 | 60° | 0.263 | 電子の 0.8 %、エネルギーの 7 %（Sect. 4.8） |
| `xu2020electron` | 1 次元 | 100 | 63（11–63 を走査） | 記載なし | 63° | 0.15（壁） | DSA。数で 7 %（Sect. 3） |
| `kumar2021nonthermal` | 1 次元 | 256 | 63 | 各 ≈1.8 | 75° | 0.15 | 両粒子種で DSA が始まる（Sect. 3） |

**斜め・準平行（θ_Bn < 60°）**

| 文献 | 次元 | m_i/m_e | M_A | β | θ_Bn | v_sh/c | 結果 |
|---|---|---|---|---|---|---|---|
| `morris2022preacceleration` | 2 次元、面内（63°、74.3° は面外も） | 64 | 21 | 1 | 30–74.3° | 0.264 | 反射率 5 %（30°）→ 0.05 %（63°）→ 0（74.3°） |
| `park2015simultaneous` | 1 次元 | 100, 400 | 20 | 記載なし | 30° | 0.1（流速） | DSA。K_ep = 3.8e-3 |
| `kato2015particle` | 1 次元 | 30 | ≈28 | 1 | 30° | 0.46 | 注入は 1 回のバースト。エネルギーの指数 ≈3.0 |
| `fang2019early` | 1 次元 | 30 | 16.5 | 記載なし | 15, 30, 45° | 0.2 | 3 角度すべてで尾。DSA は "not well developed"（Sect. 4） |
| `crumley2019kinetic` | 2 次元、面内 | 64 | 15 | 1 | 10°, 55° | 0.83 | 10°：p^(−4.2)。55°：尾なし |
| `shalaby2022mechanism` | 1 次元 | 1836, 100 | 5.3, 21.3 | ~1e-3 | 0° | 0.133 | 低い M_A と現実の質量比でだけ効率がよい（Sect. 6） |
| `arbutina2021nonlinear` | 薄い 2 次元 | 16–100 | 11–16 | 記載なし | 0° | 0.33–0.4 | p^(−4.2)。η ≈ 1e-4（Sect. IV） |
| `gupta2024electron` | 準 1 次元、16 本 | 16–1836 | 5–40 | 0.06–64 | 30° | 0.067–0.267 | 非熱的は数で ≲ 0.5 %。M_A, M_s ≳ 10 で p⁻⁴ |
| `gupta2025speed` | 準 1 次元 | 100 | 20, 5 | 記載なし | 30° | 0.133 | 閾値：速さが v_sh を超えること |
| `zekovic2025slams` | 1 次元、2 次元面内 | 32 | 20–300 | 記載なし | 0° | 0.033–0.33 | DSA の閾値より下で傾き 5–5.7 |
| `jikei2026dependence` | 2 次元、面内 | 100 | 100 付近で領域が切り替わる | 記載なし | 20° | 準相対論的 | 衝撃波エネルギーの 1.5 %（Bell）に対し 15 %（Weibel） |

**ハイブリッド場または与えた場の中のテスト電子**

| 文献 | 方法 | M_A | θ_Bn | 結果 |
|---|---|---|---|---|
| `burgess2006simulations` | 2 次元ハイブリッド、テスト電子 | 5.1, 2.1 | 80–88° | 波打ちが最大エネルギーをほぼ 10 倍にする（Fig. 5, Sect. 4） |
| `trotta2019electron` | 2 次元・3 次元ハイブリッド、テスト電子 | 2.9–6.6 | 80–87° | M = 3.5 付近に閾値（Sect. 5） |
| `guo2010effect` | 乱流を与えた 2 次元ハイブリッド、テスト電子 | 4.0 | 平均 60–90° | 90° で初期の 100 eV の 200–300 倍まで（Fig. 4） |
| `guo2015acceleration` | 3 次元テスト電子、与えた衝撃波と乱流 | ≫ 1 を仮定 | 平均 0–90° | 10 keV を超える割合は、0–90° にわたり 0–9 %（弱い乱流）から 7–11 %（強い乱流）（Table 1） |

**低マッハ数の参照例**（SNR 相当ではない）：`guo2014nonthermal`、2 次元、M_s = 3、β_p = 20、θ_Bn = 63°。
衝撃波ドリフトと、戻ってくる電子が駆動する波によって、加速の繰り返しが起き、電子の約 15 % が尾に入る
（Sect. 5, Fig. 4c）。

**「非熱的」の定義は文献によって違う。** そのため、上の割合は直接には比べられない：
- マクスウェル分布のフィットからの超過（`bohdan2017electron`）。
- 熱的運動量の 5 倍より上（`gupta2024electron` Eq. 3）。
- 熱的ピークのエネルギーの 5 倍より上（`xu2020electron`）。
- 4 m_e c より上（`kumar2021nonthermal`）。
- ピーク運動量の 3 倍より上（`jikei2026dependence`）。

## 6. θ_Bn への依存、平面でない衝撃波面、全体の幾何

### 6.1 どの角度が計算されているか
- 90° のみ：多次元の研究の大部分（5 節）。
- 一つの研究の中での走査：
  - M_A = 7、m_i/m_e = 400 で 45–90°。2 次元、ホイッスラーの領域（`riquelme2011electron` Sect. 4.1）。
  - M_A = 21 で 30–74.3°。2 次元、反射とフォアショックのみ、短い計算（`morris2022preacceleration`）。
  - 1 次元で 70–85°（`amano2007electron`）。
  - 1 次元で 15–45°（`fang2019early`）。
- それ以外は単一の角度：0°、10°、20°、30°、55°、60°、63°、74.3°、75°。
- 唯一の 3 次元計算は 74.3° で、フォアショックができるには "too short" である（`bohdan2023electron`
  Sect. 4）。
- パラメータを固定したまま、θ_Bn の全域で、熱的な集団から DSA まで電子を追った研究はない。
  45° と 63° の間の斜めの範囲を、DSA まで系統的に走査した研究もない。

### 6.2 依存性について述べられていること
- 衝撃波ドリフトは、θ_Bn が大きいほど利得が増え、反射される電子が減る（`amano2007electron` Sect. 1）。
  支配する量は M_A / cos θ_Bn である（`katou2019theory` Eq. 21、`amano2022theory` Eq. 14、
  `wang2026generation` Eq. 22）。
- M_A = 10 では依存性は単調でない：75° と 80° は限界的、85° は SSDA に必要なホイッスラー波に対して不安定、
  87° では再び弱い（`wang2026generation` Sect. IV C）。
- ホイッスラーの領域では、最も硬いスペクトルは M_A = 7 で 70–75°、M_A = 14 で 60° に現れる。90° では尾が
  できない（`riquelme2011electron` Sects. 4.1, 4.4）。電子が上流へ逃げる条件は
  100 cos²θ_Bn ≳ M_A / (m_i/m_e)^(1/2) である（Sect. 5）。
- 検証されていない推測：
  - 準平行衝撃波は両粒子種を加速し、準垂直衝撃波は電子だけを加速する（`park2015simultaneous`、
    Conclusions）。
  - 電子にとっての最適な角度は 60° と 70° の間にある（`xu2020electron` Sect. 3）。
  - 前駆領域が成長すると、輸送は初期の角度によらなくなる（`kumar2021nonthermal` Sect. 4）。
- SNR 向けのモデルでの上限：85° 付近で電子は逃げられなくなり、88° 付近で最大エネルギーが制限される
  （`amano2007electron` Sect. 3.5）。

### 6.3 陽子が駆動する波がどこにあるかは陽子で決まる
- ハイブリッド計算では、イオンの加速効率は約 40° より下で 5–15 %、45° より上で下がり、60° より上では
  ゼロになる（`caprioli2014simulations` Fig. 3, Sect. 7）。上流の δB/B₀ は準平行衝撃波で 3–6 に達し、
  3 次元の 80° では 0.1 未満にとどまる（Figs. 4, 14）。
- したがって、一様磁場中の球面では、陽子が駆動する乱流は準平行の極冠にしかない（本調査の解釈）。
- 一様磁場の場合、θ_Bn < 45° の極冠は球面の 29 %、45° から 75° の帯は 45 %、75° より上の帯は 26 % である
  （本調査の計算）。
- これには異論がある：別の手法では準垂直衝撃波でもイオンが加速される（van Marle et al. 2018、
  Reville & Bell 2013。未読。`pohl2020pic` Sect. 5.3 経由）。Orusa & Caprioli 2023 は 3 次元ハイブリッド計算で
  それを見いだしている（未読。`mcginness2026sn` が引用）。
- 陽子が注入されないままの電子注入は、63° で見られている（`xu2020electron` Sect. 3）。

### 6.4 局所的な角度：波打ち、乱流、自己生成磁場
一つの衝撃波面に沿った局所的な θ_Bn は、平均とは違う。
- **波打ち（リップル）。** 準垂直衝撃波面のイオンスケールの波打ちが、層の中で電子を散乱する。局所的な角度の
  広がりは、M = 2.9 で 1.4°、M = 6.6 で 15.2° である（`trotta2019electron` Fig. 2）。電子は、面が局所的に
  垂直な場所から入り、そうでない場所から出ることができる（Sect. 5）。2 次元の波打ちは捕捉を過大にし、
  3 次元ではスペクトルの落ち方が緩やかになる（Figs. 4, 9, 10）。これらはテスト電子を入れたハイブリッド計算で
  ある。
- **あらかじめある乱流。** 蛇行する磁力線が、電子を何度も衝撃波面に連れ戻す（`guo2010effect`）。乱流が
  強いと、平均磁場が平行な衝撃波でも、局所的な磁場が準垂直な場所で電子が加速される
  （`guo2015acceleration` Table 1）。
- **準平行衝撃波での自己生成磁場。**
  - 陽子が駆動する波が、局所的な角度を約 80° まで上げる（`kato2015particle` Sects. 3.1, 4）。
  - arccos(v_sh/c) を超える部分は超光速になる（`gupta2024electron` Sect. 4）。
  - 磁化された電子にとって、平行衝撃波は準垂直衝撃波に似てくる（`zekovic2025slams` Sects. 4.4, 5、
    `jikei2026dependence` Sect. 3.3）。
  - 同じ考えが、若い残骸の動径磁場との関係としても提案されている（`amano2020observational` p. 5）。

### 6.5 大きなスケールで曲がった衝撃波面
- SNR 相当のパラメータで、球面などの曲がった衝撃波での電子注入を扱った運動論的シミュレーションは
  見つからなかった。
- 曲がった衝撃波の文献は地球のバウショックに関するもので、読んでいない（10.2 節）：
  - Krauss-Varban & Burgess 1991：散乱のない曲がった衝撃波。90° 付近で入った電子は、わずかに低い θ_Bn で
    出ていく（`burgess2006simulations` 経由）。
  - Vandas 1995, 2001：有限の厚みを持つ曲がった衝撃波。断熱理論は観測されたスペクトルを説明できない。
  - Krauss-Varban 1994：曲がった衝撃波での散乱。
  - Savoini & Lembège 2001：曲がった衝撃波とその電子フォアショックの全粒子シミュレーション。
- バウショックでは、曲率によって、衝撃波ドリフトのエネルギーが散乱の有無によらず約 10 keV に制限される。
  これは "not an intrinsic limit" とされる（`katou2019theory` Sect. 4）。`amano2022theory` はこの上限を式で
  与えている（Eq. 22, Sect. 4.5）。
- 加速された電子は磁力線に沿って移動し、"due to the change in local θ_Bn" によって逃げる
  （`katou2019theory` Sect. 4）。バウショックで最もエネルギーの高い上流電子は、磁力線の接点から来る
  （Anderson et al. 1979、未読。`amano2022nonthermal` Sect. 3.1 経由）。
- 残骸全体は、MHD や半解析的な手法に、仮定した ε(θ_Bn) を入れてモデル化されている（3.3 節）。その場合、
  明るさのピークは磁場の軸から約 35° のところに来る。磁場の圧縮が、仮定した効率の低下より速く θ とともに
  増えるためである（`winner2020evolution` Sect. 3.2, Fig. 5）。

### 6.6 非定常性
- 衝撃波の再形成の周期は 1.5–2 Ω_ci⁻¹ である（`bohdan2017electron` Sect. 3.2、`amano2007electron`
  Sect. 3.4、`fulat2026pic` Sect. 3.1）。
- Buneman 波の強さ、サーフィン、リコネクション、反射される割合は、この周期に従って変わる
  （`bohdan2017electron` Fig. 15、`kumar2021nonthermal` Sect. 4、`gupta2025speed` App. C）。
- 注入は間欠的になりうる：`kato2015particle` では 1 回のバーストだけである（Figs. 15, 16）。

## 7. 文献が認めている数値上の制約

- **質量比。** 閾値のスケーリングが機構ごとに違う：サーフィンは (m_i/m_e)^(2/3)、ホイッスラーと SSDA は
  (m_i/m_e)^(1/2)、Buneman 不安定性は (β_e m_i/m_e)^(1/2)。そのため、計算が、実際の衝撃波とは違う領域に
  入ってしまうことがある（本調査の解釈）。
  - m_i/m_e = 16 では、電子のエネルギーの割合が 10 倍以上大きく出る（`crumley2019kinetic` App. A）。
  - `shalaby2022mechanism` の平行衝撃波の機構は、100 では消える。
  - `gupta2024electron` は、割合が質量比によらないとする。
- **衝撃波速度。** PIC 計算の多くは 0.1 c を超える速度を使い、現実の 10 倍である（`pohl2020pic`
  Sects. 2.2, 5）。1 次元計算の一部は 0.03–0.07 c まで下げている（`zekovic2025slams`、`gupta2024electron`）。
  その結果：
  - θ_Bn の超光速の範囲が広くなる（2 節）。
  - サーフィンの反射 1 回がすでに相対論的になる。これは 0.01 c では成り立たない（`amano2022nonthermal`
    Sect. 4.5）。
  - 注入までの隔たりが、10⁵–10⁶ 倍ではなく 10³ 倍になる（`bohdan2023electron` Sect. 4.3）。
  - 衝撃波速度に対するスケーリングは "not entirely clear" である（`bohdan2019kineticb` Sect. 4）。
- **次元。**
  - 1 次元には平面の Buneman 波があり、Weibel モードはない。
  - 磁場が面外の 2 次元には Buneman モードがあり、ホイッスラー波と波打ちは排除される。
  - 磁場が面内の 2 次元にはホイッスラー波、Weibel モード、リコネクションがあり、Buneman モードは
    抑えられる。
  - すべてを含むのは 3 次元だけである（`riquelme2011electron` Appendix、`matsumoto2017electron`、
    `bohdan2017electron` Sect. 4）。
  - "Only 3D simulations can provide a complete picture"（`pohl2020pic` Sect. 5.3）。
- **計算時間と箱。** 2 次元では数から数十 Ω_ci⁻¹。フォアショックには 125–300 が必要である（4.2 節）。
  いくつかの箱は狭すぎて、波打ちが入らない（`bohdan2017electron` Sect. 3.2、`arbutina2021nonlinear`
  Sect. II）。
- **数値的な抵抗。** PIC では速い粒子に数値的な抵抗が働く（Kato 2013、arXiv:1312.5507、未読。
  `kato2015particle`、`xu2020electron`、`kumar2021nonthermal` が引用）。

## 8. 食い違い、未解決の問題、抜けている点

### 8.1 文献の間の食い違い
1. **現実のパラメータでサーフィンは効くか。** マッハ数の閾値より上では効く（`matsumoto2012electron`、
   `matsumoto2017electron`）。効かない（`ohira2007absence`、`riquelme2011electron`、`bohdan2019kineticb`）。
2. **SNR 衝撃波にはどのマッハ数が必要か。** ホイッスラーには M_A ≲ 20（`riquelme2011electron`）。
   サーフィンには M_A ≳ (m_i/m_e)^(2/3)（`matsumoto2012electron`）。
3. **高マッハ数での衝撃波ドリフト。** "cannot operate efficiently"（`guo2014nonthermal` Sect. 5）に対し、
   M_A / cos θ_Bn が大きいほど有利（`katou2019theory`、`amano2020observational`）。
4. **平行衝撃波。** 陽子が駆動する乱流と反射（`park2015simultaneous`、`gupta2024electron`）。現実の質量比と
   低い M_A を必要とするイオンサイクロトロン波（`shalaby2022mechanism`）。高 M_A での急峻化した構造
   （`zekovic2025slams`）。
5. **小さくした質量比。** 信頼できる（`gupta2024electron`、`arbutina2021nonlinear`）。信頼できない
   （`shalaby2022mechanism`）。
6. **60° より上でのイオン。** 63° では急な尾だけ（`xu2020electron`）。75° ではイオンの尾ができる
   （`kumar2021nonthermal`）。
7. **電子にとって準平行か準垂直か。** PIC での効率は 30° より 63° のほうが高い（`xu2020electron` Sect. 3、
   `gupta2024electron` Sect. 7）。SN 1006 のモデルは逆を必要とする（`winner2020evolution`）。電波残骸の集団は
   準垂直を支持する（`west2017connection`）。
8. **注入の閾値**そのもの（2 節）。
9. **バウショックでエネルギーを制限するもの：** 曲率（`katou2019theory`）か、散乱
   （`amano2020observational`）か。

### 8.2 文献が挙げている未解決の問題
- 電子の効率の、角度に対する関数形（`winner2020evolution` Sect. 2.2）。
- サーフィン、SSDA、リコネクションの相対的な役割。散乱率。SSDA が Weibel 乱流の中で生き残るか
  （`amano2022nonthermal` Sects. 5.4, 6）。
- 衝撃波速度、質量比、マッハ数に対するスケーリング。多次元の平行衝撃波（`bohdan2023electron` Sect. 6）。
- 注入効率が衝撃波の条件にどう依存するか。ホイッスラー波で加速された電子から p⁻⁴ への移行
  （`marcowith2016microphysics` Sect. 3.4.1）。
- 輸送がより複雑になる、より斜めの衝撃波（`gupta2025speed` Sect. IV）。
- 準垂直衝撃波での電子の DSA の、多次元での確認（`xu2020electron` Sect. 4）。
- 完全に運動論的な斜め衝撃波の計算。3 次元（`pohl2020pic` Sect. 5.3）。
- 上流の乱流のどの部分が効くのか：磁場による散乱か、衝撃波面の変形か（`fulat2026pic` Sect. 3.3）。

### 8.3 このテーマから見て抜けている点（本調査の解釈）
- **曲がった面がない。** 大きなスケールで曲がり、θ_Bn が面に沿って変わる衝撃波面での電子注入は、見つけた
  文献の中では、SNR 相当のパラメータでの運動論的シミュレーションがない。6.5 節のバウショックの文献が直接の
  先行研究であり、読んでいない。それを読むまで、曲率についての新規性の主張は確実ではない。
- **一組のパラメータでの走査がない。** 注入の θ_Bn 依存は、質量比、マッハ数、次元、定義の違う研究から、
  断片的に知られているだけである。
- **全体モデルはまさにその関数を必要としている。** 電子の ε(θ_Bn) は 3.3 節のどのモデルでも仮定であり、
  SN 1006 の観測は、仮定した形のどれを選ぶかに使われている。
- **面に沿った輸送。** 局所的な θ_Bn の変化によって電子が元の領域を離れる、という記述は、バウショックに
  ついてしかない（6.5 節）。
- **領域。** SNR 衝撃波は M_A が数十から数百で、上のシミュレーションでは Weibel 媒介である。M_A が数程度の
  衝撃波は、`trotta2019electron`、`burgess2006simulations`、`guo2014nonthermal` の領域にある。球面の計算が
  どの領域に入るかで、どの文献と比べられるかが決まる。
- **超光速の帯。** シミュレーションの速度では、球面のかなりの部分が超光速になる（2 節）。SNR の速度では
  2 % 未満である。

これらは、テーマの研究の問いを決めるものではない。その判断はオーナーのものである。

## 9. 検索の記録と、この調査の限界

**検索**（2026-11-12、arXiv API、要旨と題名。各検索で新しい順に 60 件。既知でない結果は 238 件）。
検索語：
- "electron injection" と shock
- "electron acceleration"、"supernova remnant"、"particle-in-cell"
- "quasi-parallel"、electron、shock、"particle-in-cell"
- "shock drift acceleration" と electron
- "shock surfing" と electron
- "curved shock"、"shock curvature"、"spherical shock" のいずれかと、electron、acceleration
- obliquity、electron、shock、acceleration
- "electron-to-proton"、shock、acceleration
- "SN 1006" と、"quasi-parallel"、obliquity、"polar cap" のいずれか
- "electron foreshock"
- "pre-acceleration"、electron、shock
- Weibel、electron、"supernova remnant"、shock

2 回目は著者と題名で検索した（上の文献の研究グループ、残骸の X 線偏波、双極型残骸、実験室の衝撃波、
PIC での爆発殻）。読んだ文献の参考文献欄からも文献を拾った。

**選定。** この調査のために 54 本を読んだ（10.1 節）。SNR 相当の衝撃波での電子注入と、角度依存に直接関わる
ものを選んだ。それ以前からあるノートが 5 本ある。

**読んだ深さ。** 各ノートに書いてある。総説は指定した節だけを読んだ。直接の先行研究としなかった論文は、
導入、設定、結論を読み、結果は数値を拾うために流し読みした。

**確認。** ノートを書いたあと、20 本の主要な数値を論文本文と照合した。図から数値を読み取ることは
していない。ノートがそう断っている場合（`amano2007electron` Fig. 7）は除く。

**限界。**
- 到達できたのは arXiv だけである。そこにない文献は読んでいない：1980 年代・1990 年代の理論、
  J. Geophys. Res. のバウショックの文献、Matsumoto et al. 2015（Science）、いくつかの観測論文。これらは
  10.2 節にある。
- 雑誌の巻とページは、arXiv のメタデータか、ここで読んだ他の文献の参考文献欄による。どちらにもない場合は
  「未確認」と書いた。
- 曲がった衝撃波の検索では、SNR 相当の衝撃波についての運動論的な研究は見つからなかった。1 回の検索は、
  存在しないことの証明ではない。
- 扱っていないもの：`guo2014nonthermal` 以外の銀河団の低マッハ数衝撃波、相対論的衝撃波、実験室実験、
  電子加熱だけの研究、惑星間衝撃波のその場観測。
- 被引用（主要な文献を誰が引用しているか）は、キーワードでしか調べていない。引用索引は使っていない。

## 10. 参考文献

### 10.1 読んだ文献（ノートは `research/literature/`）

観測と、残骸全体のモデル：
- `vink2012supernova` — J. Vink 2012, Astron. Astrophys. Rev. 20, 49（論文番号は、他の文献の参考文献欄の記載による）. arXiv:1112.0576. DOI 10.1007/s00159-011-0049-1.
- `reynolds2004microphysics` — S. P. Reynolds 2004, Adv. Space Res.（COSPAR 会議録。巻は未確認）.
  arXiv:astro-ph/0308217.
- `morlino2012strong` — G. Morlino & D. Caprioli 2012, A&A 538, A81. arXiv:1105.6342.
  DOI 10.1051/0004-6361/201117855.
- `bocchino2011constraints` — F. Bocchino, S. Orlando, M. Miceli & O. Petruk 2011, A&A 531, A129.
  arXiv:1105.2689. DOI 10.1051/0004-6361/201016341.
- `reynoso2013radio` — E. M. Reynoso, J. P. Hughes & D. A. Moffett 2013, AJ 145, 104. arXiv:1302.4678.
  DOI 10.1088/0004-6256/145/4/104.
- `west2017connection` — J. L. West, S. Safi-Harb & G. Ferrand 2017, A&A 597, A121. arXiv:1609.02887.
  DOI 10.1051/0004-6361/201628079.
- `zhou2023magnetic` — P. Zhou et al. 2023, ApJ（受理済み。巻は未確認）. arXiv:2309.01879.
- `giuffrida2022supernova` — R. Giuffrida et al. 2022, Nat. Commun. 13, 5098. arXiv:2208.14491.
- `mcginness2026sn` — E. McGinness, R. Diesing, D. Caprioli & F. Acero 2026, ApJ（巻は未確認）.
  arXiv:2606.12111. DOI 10.3847/1538-4357/aea3f3.
- `winner2020evolution` — G. Winner, C. Pfrommer, P. Girichidis, M. Werhahn & M. Pais 2020, MNRAS 499, 2785.
  arXiv:2006.06683. DOI 10.1093/mnras/staa2989.

準垂直衝撃波の PIC：
- `amano2007electron` — T. Amano & M. Hoshino 2007, ApJ 661, 190. arXiv:astro-ph/0612204. DOI 10.1086/513599.
- `amano2009electron` — T. Amano & M. Hoshino 2009, ApJ 690, 244. arXiv:0805.1098.
  DOI 10.1088/0004-637X/690/1/244.
- `ohira2007absence` — Y. Ohira & F. Takahara 2007, ApJ 661, L171. arXiv:0705.2061.
- `matsumoto2012electron` — Y. Matsumoto, T. Amano & M. Hoshino 2012, ApJ 755, 109. arXiv:1204.6312.
  DOI 10.1088/0004-637X/755/2/109.
- `matsumoto2013electron` — Y. Matsumoto, T. Amano & M. Hoshino 2013, Phys. Rev. Lett. 111, 215003.
  arXiv:1310.8367. DOI 10.1103/PhysRevLett.111.215003.
- `matsumoto2017electron` — Y. Matsumoto, T. Amano, T. N. Kato & M. Hoshino 2017, Phys. Rev. Lett. 119, 105101.
  arXiv:1709.03673. DOI 10.1103/PhysRevLett.119.105101.
- `kato2010nonrelativistic` — T. N. Kato & H. Takabe 2010, ApJ 721, 828. arXiv:1008.0265.
  DOI 10.1088/0004-637X/721/1/828.
- `riquelme2011electron` — M. A. Riquelme & A. Spitkovsky 2011, ApJ 733, 63. arXiv:1009.3319.
  DOI 10.1088/0004-637X/733/1/63.
- `wieland2016nonrelativistic` — V. Wieland, M. Pohl, J. Niemiec, I. Rafighi & K.-I. Nishikawa 2016, ApJ 820, 62.
  arXiv:1602.05064. DOI 10.3847/0004-637X/820/1/62.
- `bohdan2017electron` — A. Bohdan, J. Niemiec, O. Kobzar & M. Pohl 2017, ApJ 847, 71. arXiv:1708.05528.
  DOI 10.3847/1538-4357/aa872a.
- `bohdan2019kinetic` — A. Bohdan et al. 2019, ApJ 878, 5. arXiv:1904.13153. DOI 10.3847/1538-4357/ab1b6d.
- `bohdan2019kineticb` — A. Bohdan et al. 2019, ApJ 885, 10. arXiv:1909.05294. DOI 10.3847/1538-4357/ab43cf.
- `bohdan2020kinetic` — A. Bohdan et al. 2020, ApJ 893, 6. arXiv:2003.01879. DOI 10.3847/1538-4357/ab7cd6.
- `bohdan2021magnetic` — A. Bohdan et al. 2021, Phys. Rev. Lett. 126, 095101. arXiv:2102.04328.
  DOI 10.1103/PhysRevLett.126.095101.

斜め衝撃波、フォアショック、乱流のある上流：
- `morris2022preacceleration` — P. J. Morris, A. Bohdan, M. S. Weidl & M. Pohl 2022, ApJ 931, 129.
  arXiv:2204.11569. DOI 10.3847/1538-4357/ac69c7.
- `morris2023preacceleration` — P. J. Morris et al. 2023, ApJ 944, 13. arXiv:2301.00872.
  DOI 10.3847/1538-4357/acaec8.
- `fulat2023kinetic` — K. Fulat, A. Bohdan, G. Torralba Paz & M. Pohl 2023, ApJ 959, 119. arXiv:2310.14388.
- `fulat2026pic` — K. Fulat et al. 2026, ApJ に受理済み（巻は未確認）. arXiv:2605.21986.

準平行衝撃波と注入の閾値：
- `park2015simultaneous` — J. Park, D. Caprioli & A. Spitkovsky 2015, Phys. Rev. Lett. 114, 085003.
  arXiv:1412.0672. DOI 10.1103/PhysRevLett.114.085003.
- `kato2015particle` — T. N. Kato 2015, ApJ 802, 115. arXiv:1407.1971.
- `crumley2019kinetic` — P. Crumley, D. Caprioli, S. Markoff & A. Spitkovsky 2019, MNRAS 485, 5105.
  arXiv:1809.10809. DOI 10.1093/mnras/stz232.
- `xu2020electron` — R. Xu, A. Spitkovsky & D. Caprioli 2020, ApJ Lett. 897, L41. arXiv:1908.07890.
  DOI 10.3847/2041-8213/aba11e.
- `kumar2021nonthermal` — N. Kumar & B. Reville 2021, ApJ Lett. 921, L14. arXiv:2110.09939.
  DOI 10.3847/2041-8213/ac30e0.
- `fang2019early` — J. Fang, C.-Y. Lu, J.-W. Yan & H. Yu 2019, Res. Astron. Astrophys. 19, 182. arXiv:1908.08170.
- `shalaby2022mechanism` — M. Shalaby, R. Lemmerz, T. Thomas & C. Pfrommer 2022, ApJ 932, 86. arXiv:2202.05288.
  DOI 10.3847/1538-4357/ac6ce7.
- `arbutina2021nonlinear` — B. Arbutina & V. Zeković 2021, Astropart. Phys. 127, 102546. arXiv:2012.15117.
  DOI 10.1016/j.astropartphys.2020.102546.
- `gupta2024electron` — S. Gupta, D. Caprioli & A. Spitkovsky 2024, ApJ 976, 10. arXiv:2408.16071.
  DOI 10.3847/1538-4357/ad7c4c.
- `gupta2025speed` — S. Gupta, D. Caprioli & A. Spitkovsky 2025, プレプリント arXiv:2506.09134. おそらく
  ApJ Lett. 994, L34 (2025) として出版されている：`mcginness2026sn` が K_ep の式の出典としてその項目を引用している。
  同一の論文かどうかは確認していない。
- `zekovic2025slams` — V. Zeković, A. Spitkovsky & Z. Hemler 2025, ApJ 988, 40. arXiv:2408.02084.
  DOI 10.3847/1538-4357/ade05f.
- `jikei2026dependence` — T. Jikei, D. Grošelj & L. Sironi 2026, プレプリント arXiv:2607.05778.

衝撃波ドリフトの理論、太陽圏での検証、幾何：
- `katou2019theory` — T. Katou & T. Amano 2019, ApJ 874, 119. arXiv:1903.02277. DOI 10.3847/1538-4357/ab0d8a.
- `amano2020observational` — T. Amano et al. 2020, Phys. Rev. Lett. 124, 065101. arXiv:2002.06787.
  DOI 10.1103/PhysRevLett.124.065101.
- `wang2026generation` — R. Wang & T. Amano 2026, Phys. Plasmas に受理済み（巻は未確認）.
  arXiv:2603.19018.
- `raptis2025revealing` — S. Raptis et al. 2025, Nat. Commun. 16, 488. arXiv:2502.10643.
  DOI 10.1038/s41467-024-55641-9.
- `guo2014nonthermal` — X. Guo, L. Sironi & R. Narayan 2014, ApJ 794, 153. arXiv:1406.5190.
  DOI 10.1088/0004-637X/794/2/153.
- `burgess2006simulations` — D. Burgess 2006, ApJ 653, 316. arXiv:astro-ph/0610714. DOI 10.1086/508805.
- `guo2010effect` — F. Guo & J. Giacalone 2010, ApJ 715, 406. arXiv:1003.5946.
  DOI 10.1088/0004-637X/715/1/406.
- `guo2015acceleration` — F. Guo & J. Giacalone 2015, ApJ 802, 97. arXiv:1409.5854.
  DOI 10.1088/0004-637X/802/2/97.
- `trotta2019electron` — D. Trotta & D. Burgess 2019, MNRAS 482, 1154. arXiv:1808.00812.
  DOI 10.1093/mnras/sty2756.
- `caprioli2014simulations` — D. Caprioli & A. Spitkovsky 2014, ApJ 783, 91. arXiv:1310.2943.
  DOI 10.1088/0004-637X/783/2/91.

総説：
- `bohdan2023electron` — A. Bohdan 2023, Plasma Phys. Control. Fusion 65（論文番号は未確認）.
  arXiv:2211.13992. DOI 10.1088/1361-6587/aca5b2.
- `amano2022nonthermal` — T. Amano et al. 2022, Rev. Mod. Plasma Phys. 6, 29. arXiv:2209.03521.
- `marcowith2016microphysics` — A. Marcowith et al. 2016, Rep. Prog. Phys. 79, 046901. arXiv:1604.00318.
  DOI 10.1088/0034-4885/79/4/046901.
- `pohl2020pic` — M. Pohl, M. Hoshino & J. Niemiec 2020, Prog. Part. Nucl. Phys., article 103751（巻は未確認）.
  arXiv:1912.02673. DOI 10.1016/j.ppnp.2019.103751.

この調査の前からあったノート：
- `hoshino2002nonthermal` — M. Hoshino & N. Shimada 2002, ApJ 572, 880. arXiv:astro-ph/0203073.
- `amano2010critical` — T. Amano & M. Hoshino 2010, Phys. Rev. Lett. 104, 181102. arXiv:1004.3838.
- `amano2022theory` — T. Amano & M. Hoshino 2022, ApJ 927, 132. arXiv:2201.11416.

### 10.2 本文で引用したが読んでいない文献
各項目は、角括弧内の文献の参考文献欄に記載されたとおりに転記した。

残骸の観測とモデル：
- Fulbright, M. S. & Reynolds, S. P. 1990, ApJ, 357, 591 [`bocchino2011constraints`, `west2017connection`].
- Leckband, J. A., Spangler, S. R., & Cairns, I. H. 1989, ApJ, 338, 963 [`west2017connection`].
- Rothenflug, R., Ballet, J., Dubner, G., et al. 2004, A&A, 425, 121 [`bocchino2011constraints`].
- Koyama, K., et al. 1995, Nature, 378, 255 [`vink2012supernova`].
- Jokipii, J. R. 1987, ApJ, 313, 842 [`reynolds2004microphysics`].
- Cassam-Chenaï, G., Hughes, J. P., Reynoso, E. M., Badenes, C., & Moffett, D. 2008, ApJ, 680, 1180
  [`reynoso2013radio`].
- Gaensler, B. M. 1998, ApJ, 493, 781 [`bocchino2011constraints`]. arXiv にある：astro-ph/9707149.
- Völk, H. J., Berezhko, E. G., & Ksenofontov, L. T. 2003, A&A, 409, 563 [`bocchino2011constraints`,
  `winner2020evolution`]. arXiv にある：astro-ph/0306016.
- Ellison, D. C., Baring, M. G., & Jones, F. C. 1995, ApJ, 453, 873 [`bocchino2011constraints`].

サーフィン、Buneman 不安定性、Weibel 乱流：
- Buneman, O. 1958, Phys. Rev. Lett., 1, 8 [`bohdan2017electron`].
- Papadopoulos, K. 1988, Ap&SS, 144, 535 [`amano2007electron`].
- Shimada, N., & Hoshino, M. 2000, ApJ, 543, L67 [`amano2007electron`, `wieland2016nonrelativistic`].
- McClements, K. G., Dieckmann, M. E., Ynnerman, A., Chapman, S. C., & Dendy, R. O. 2001, Phys. Rev. Lett.,
  87, 255002 [`amano2007electron`].
- Matsumoto, Y., Amano, T., Kato, T. N., & Hoshino, M. 2015, Science, 347, 974 [`bohdan2020kinetic`].
- Bohdan, A., Weidl, M. S., Morris, P. J., & Pohl, M. 2022, Physics of Plasmas, 29, 052301
  [`morris2023preacceleration`, `fulat2026pic`].

衝撃波ドリフト、ホイッスラー波、太陽圏の観測：
- Wu, C. S. 1984, J. Geophys. Res., 89, 8857 [`katou2019theory`, `amano2007electron`].
- Leroy, M. M., & Mangeney, A. 1984, Annales Geophysicae, 2, 449 [`katou2019theory`].
- Krauss-Varban, D., & Wu, C. S. 1989, J. Geophys. Res., 94, 15367 [`katou2019theory`, `guo2014nonthermal`].
- Ball, L., & Melrose, D. B. 2001, PASA, 18, 361 [`park2015simultaneous`, `xu2020electron`].
- Levinson, A. 1992, ApJ, 401, 73 [`katou2019theory`, `riquelme2011electron`].
- Krasnoselskikh, V. V., Lembège, B., Savoini, P., & Lobzin, V. V. 2002, Physics of Plasma, 9, 1192
  [`riquelme2011electron`].
- Oka, M., Terasawa, T., Seki, Y., et al. 2006, Geophys. Res. Lett., 33, L24104 [`katou2019theory`].
- Oka, M., Wilson, III, L. B., Phan, T. D., et al. 2017, ApJL, 842, L11 [`katou2019theory`].
- Masters, A., Stawarz, L., Fujimoto, M., et al. 2013, Nature Physics, 9, 164 [`katou2019theory`,
  `guo2015acceleration`].
- Anderson, K. A., Lin, R. P., Martel, F., Lin, C. S., Parks, G. K., & Rème, H. 1979, Geophys. Res. Lett., 6,
  401 [`burgess2006simulations`, `amano2022nonthermal`].

曲がった衝撃波、波打つ衝撃波：
- Krauss-Varban, D., & Burgess, D. 1991, J. Geophys. Res., 96, 143 [`burgess2006simulations`,
  `katou2019theory`].
- Krauss-Varban, D. 1994, J. Geophys. Res., 99, 2537 [`burgess2006simulations`].
- Vandas, M. 1995, J. Geophys. Res., 100, 23499; Vandas, M. 2001, J. Geophys. Res., 106, 1859
  [`burgess2006simulations`, `katou2019theory`].
- Savoini, P., & Lembège, B. 2001, J. Geophys. Res., 106, 12975 [`burgess2006simulations`].
- Lembège, B., & Savoini, P. 2002, J. Geophys. Res., 107, 1037 [`burgess2006simulations`,
  `trotta2019electron`].
- Lowe, R. E., & Burgess, D. 2003, Ann. Geophys., 21, 671 [`burgess2006simulations`, `trotta2019electron`].
- Jokipii, J. R., & Giacalone, J. 2007, ApJ, 660, 336 [`guo2010effect`].

陽子と準平行衝撃波：
- Bell, A. R. 2004, MNRAS, 353, 550 [`park2015simultaneous`].
- Caprioli, D., Pop, A.-R., & Spitkovsky, A. 2015, ApJ, 798, L28 [`gupta2024electron`]. arXiv にある：1409.8291.
- Blasi, P., Gabici, S., & Vannoni, G. 2005, MNRAS, 361, 907 [`arbutina2021nonlinear`, `morlino2012strong`].
- Orusa, L., & Caprioli, D. 2023, PhRvL, 131, 095201 [`mcginness2026sn`]. arXiv にある：2305.10511.
- van Marle et al. 2018 と Reville & Bell 2013：`pohl2020pic` Sect. 5.3 に名前がある。書誌情報は転記していない。

### 10.3 検索で見つけたが読んでいない文献（今後の候補）
題名と番号は arXiv のメタデータによる。
- 2606.19979, R. Ferrazzoli et al. 2026, "Supernova Remnants in the IXPE era: a review".
- 2505.09712, P. Zhou et al. 2025, "X-ray polarization in SN 1006 southwest shows spatial variations and
  differences with the radio band".
- 1702.02054, S. Katsuda 2017, "Supernova of 1006 (G327.6+14.6)".
- 1502.03133, E. M. Schneiter et al. 2015, "3D MHD simulation of polarized emission in SN 1006".
- 0811.2319, O. Petruk et al. 2009, "Aspect angle for interstellar magnetic field in SN 1006".
- 0704.0890, S. Orlando et al. 2007, "On the Origin of Asymmetries in Bilateral Supernova Remnants".
- 1510.08536, J. L. West et al. 2016, "The Connection between Supernova Remnants and the Galactic Magnetic
  Field: A Global Radio Study of the Axisymmetric Sample".
- 1805.00128, M. Pais et al. 2018, "The effect of cosmic-ray acceleration on supernova blast wave dynamics".
- 1905.07414, R. Diesing & D. Caprioli 2019, "On the Spectrum of Electrons Accelerated in Supernova Remnants".
- 2103.02375, P. Cristofari et al. 2021, "Cosmic ray protons and electrons from supernova remnants".
- 1209.1466, K. Maeda 2012, "Injection and Acceleration of Electrons at A Strong Shock: Radio and X-ray Study
  of Young Supernova 2011dh".
- 0908.3217, G. Morlino 2009, "Cosmic-ray electron injection from the ionization of nuclei".
- 1011.4775, V. N. Zirakashvili & F. A. Aharonian 2010, "Radioactivity and Electron Acceleration in Supernova
  Remnants".
- astro-ph/0002347, M. E. Dieckmann et al. 2000, "Electron acceleration due to high frequency instabilities at
  supernova remnant shocks".
- 0808.3195, Y. Ohira & F. Takahara 2008, "Oblique Ion Two-Stream Instability in the Foot Region of a
  Collisionless Shock".
- 1004.2771, N. Shimada, M. Hoshino & T. Amano 2010, "Structure of a strong supernova shock wave and rapid
  electron acceleration confined in its transition region".
- 2008.05920, A. Bohdan et al. 2020, "Kinetic simulation of nonrelativistic perpendicular shocks of young
  supernova remnants. IV. Electron heating".
- 2308.15835, G. Torralba Paz et al. 2023, "Prediction and Anomaly Detection of accelerated particles in PIC
  simulations using neural networks".
- 1209.3938, J. Niemiec et al. 2012, "Nonrelativistic parallel shocks in unmagnetized and weakly magnetized
  plasmas".
- 1409.7393, X. Guo, L. Sironi & R. Narayan 2014, "Non-Thermal Electron Acceleration in Low Mach Number
  Collisionless Shocks. II.".
- 2202.02689, J. Fang et al. 2022, "Kinetic Simulation on Electron, Proton and Helium Acceleration in a
  Nonrelativistic Quasiparallel shock".
- 2312.13365, S. Gupta, D. Caprioli & A. Spitkovsky 2023, "Return Currents in Collisionless Shocks".
- 2312.07933, T. Jikei et al. 2023, "Enhanced Magnetic Field Amplification by Ion-Beam Weibel Instability in
  Weakly Magnetized Astrophysical Shocks".
- 2512.03169, T. Jikei et al. 2025, "Magnetic Field Amplification and Particle Acceleration in Weakly Magnetized
  Trans-relativistic Electron-ion Shocks".
- 2507.13436, L. Orusa et al. 2025, "The role of three-dimensional effects on ion injection and acceleration in
  perpendicular shocks".
- 2205.01715, J. Ng et al. 2022, "Electron-scale reconnection in three-dimensional shock turbulence".
- 2001.04945, K. V. Lezhnin et al. 2020, "Kinetic simulations of electron pre-energization by magnetized
  collisionless shocks in expanding laboratory plasmas".
- 2107.00508, O. Kobzar et al. 2021, "Electron Acceleration at Rippled Low-Mach-number Shocks in High-beta
  Collisionless Cosmic Plasmas".
- 1602.08170, X. Kong et al. 2016, "Electron Acceleration at a Coronal Shock Propagating Through a Large-scale
  Streamer-like Magnetic Field".
- 2112.12114, B. D. Keenan et al. 2021, "Hybrid Particle-in-Cell Simulations of Electromagnetic Coupling and
  Waves From Streaming Burst Debris" （エンジンの blast-wave 例の先行研究）.
