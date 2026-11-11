[English](Home) · **日本語**

# pypic の example: 厳密解・発表済みの結果と照合したシミュレーション

[pypic](https://github.com/ytworks/pypic) は、GPU で高速化した、複数 GPU 対応の 3 次元特殊相対論的な電磁粒子
（PIC）コードです。陽解法（ブネマン–ボリスのプッシュ、Esirkepov の電荷保存型割り付け、Yee 格子の場）と、
全エネルギーを厳密に保存する陰解法を備えています。すべての run でガウスの法則と ∇·B = 0 が丸め誤差の範囲で
成り立ちます。

以下の各 example は、解析解またはよく知られた発表済みの結果があるプラズマの問題です。各ページには、実験条件、
物理現象のアニメーション、シミュレーションと解を比べる主要なチャートを載せています。数値・図・アニメーションは
すべて run 自身の出力から
[`scripts/build_wiki.py`](https://github.com/ytworks/pypic/blob/main/scripts/build_wiki.py) で生成しており、
手で写した値はありません。

| example | 物理 | 比較対象 | 結果 |
|---|---|---|---|
| [二流体不安定性とブネマン不安定性](Two-Stream-Instability-ja) | 静電的なビーム不安定性（e-i、e-p、e-i-p、陰解法、ブネマン） | 暖かい流体の線形分散関係 | 成長率 {{ two_stream:r.buneman.relative_error | +.1% }} … {{ two_stream:r.ei_implicit.relative_error | +.1% }}、陰解法のエネルギー誤差 {{ two_stream:r.ei_implicit.energy_error_max | .0e }} |
| [熱的プラズマ波](Plasma-Waves-ja) | 熱雑音の ω–k スペクトル: ラングミュア波、光、O/X モード、バーンスタイン波、ホイッスラー、R/L 波、イオン音波 | 運動論的・冷たい磁化プラズマの分散関係 | 誤差の中央値 {{ plasma_waves:r.perpendicular.branches.Ez: O mode.median_rel_error | .2% }}（O モード）… {{ plasma_waves:r.unmagnetized.branches.Ex: Langmuir (kinetic).median_rel_error | .1% }}（ラングミュア波） |
| [ラングミュア振動](Landau-Damping-ja) | ランダウ減衰、非線形捕捉、Dawson 振動、砕波、相対論的振動 | 運動論的な根、冷たいシートの厳密解 | 減衰率 {{ landau_damping:r.landau.gamma_rel_error | +.1% }}、砕波時刻 {{ landau_damping:r.wavebreaking.rel_error | +.1% }} |
| [電流フィラメンテーション](Filamentation-Instability-ja) | 対向ビームのワイベル不安定性 | 冷たい相対論的分散関係 | 成長率 {{ filamentation:r.pair.median_rel_error_m2_10 | .1% }}（モードについての中央値） |
| [テスト粒子軌道](Particle-Orbits-ja) | ジャイロ運動、E×B、双曲線運動、E > B、平面波中の 8 の字 | 相対論的な厳密軌道 | 誤差 ≤ {{ particle_orbits:r.gyration.max_rel_position_error | .2% }}、2 次精度（比 {{ particle_orbits:r.plane_wave.convergence_ratio | .2f }}〜{{ particle_orbits:r.gyration.convergence_ratio | .2f }}） |
| [真空へのプラズマ膨張](Plasma-Expansion-ja) | 膨張シースによるイオン加速（TNSA） | Mora（2003）の自己相似解 | 速度分布 {{ plasma_expansion:r.isothermal_hydrogen.velocity_profile_median_rel_error | .1% }}、E_max は {{ plasma_expansion:r.isothermal_hydrogen.E_max_median_rel_error | .0% }} 低い（運動論的電子） |
| [相対論的電子陽電子衝撃波](Relativistic-Shock-ja) | γ₀ = 15 の無衝突衝撃波、フェルミ加速 | 相対論的ジャンプ条件（Synge 気体、2 次元・3 次元） | β_sh {{ relativistic_shock:r.pair_3d.beta_sh_rel_error | +.1% }}、n₂/n₁ {{ relativistic_shock:r.pair_3d.compression_rel_error | +.1% }}（3 次元） |
| [相対論的磁気リコネクション](Magnetic-Reconnection-ja) | σ = 10 のペアプラズマのリコネクション、プラズモイド、べき乗則の粒子 | Harris 平衡、Werner et al.（2016） | p = {{ reconnection:r.pair_sigma10_large.power_law_index | .2f }}（Werner 2016 の値は {{ reconnection:r.pair_sigma10_large.werner2016_index | .2f }}） |
| [磁化惑星のバウショック](Bow-Shock-ja) | 磁化した風の中の双極子の磁気圏界面とバウショック（開いた 3 次元の箱、吸収球） | 圧力平衡、垂直 MHD 衝撃波と Spreiter et al.（1966）の離れ距離 | 磁気圏界面 {{ bow_shock:r.gate.magnetopause | .2f }} R_mp、衝撃波 {{ bow_shock:r.gate.shock | .2f }}（参照値 {{ bow_shock:r.gate.shock_spreiter | .2f }}）R_mp |
| [磁化天体のまわりの爆風](Blast-Wave-ja) | 双極子をもつ吸収球のまわりで、磁化プラズマの中へ膨張するデブリの殻（開いた 3 次元の箱） | 速い磁気音波の速度、Le et al.（2021）のデブリの広がりのスケーリング | 波面は速い磁気音波の速度の {{ blast_wave:r.gate.front_mach_fast | .2f }} 倍、デブリの広がりはスケーリングの {{ blast_wave:r.gate.debris_extent_ratio | .2f }} 倍 |

## ギャラリー

| | | |
|---|---|---|
| <img src="images/two_stream/ei.gif" width="280"><br>二流体: 位相空間の渦 | <img src="images/filamentation/pair.gif" width="280"><br>フィラメンテーション: 電流フィラメント | <img src="images/reconnection/pair_sigma10_large_run2.gif" width="280"><br>リコネクション: プラズモイド |
| <img src="images/landau_damping/wavebreaking.gif" width="280"><br>砕波: シートの交差 | <img src="images/relativistic_shock/pair_2d.gif" width="280"><br>相対論的衝撃波 | <img src="images/particle_orbits/orbits.gif" width="280"><br>テスト粒子軌道 |
| <img src="images/bow_shock/gate.gif" width="280"><br>バウショック: イオン密度 | <img src="images/blast_wave/gate.gif" width="280"><br>爆風: 磁場、周囲とデブリのイオン | |

## ページの作り方

- ゲート規模の設定は、クリーンなコミットから 8 × RTX PRO 6000 Blackwell GPU で新たに実行しています。長時間の
  本番実行（衝撃波、リコネクション、膨張）は `bench/results/` に記録済みの run です。各ページには、使った run の
  コミット・run id・GPU を載せています。
- `python scripts/build_wiki.py run` で未実行の設定を実行し、`collect` で全 run を解析し（図、GIF、数値を
  `docs/wiki/data.json` へ）、`render` で `docs/wiki/` のページテンプレートを埋めます。
- 各手順に必要なもの: `render` はチェックイン済みの `docs/wiki/data.json` とテンプレートだけで動きます。`run` が
  作るのはゲート規模の新規 run だけです。`collect` は記録済みの本番 run も読みますが、その出力ディレクトリ（各
  ページの「データの出所」の表に記載）は git に入っていないので、計算したマシンから `output/` にコピーして
  ください。`collect` は記録済みの各ディレクトリを evidence の記録（run id、またはコミットと計算時間）と照合し、
  欠けている場合や別の run の場合は停止します。
- 単位はすべて規格化単位です: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、運動量は u = γv/c、温度は θ = T/(mc²)。
