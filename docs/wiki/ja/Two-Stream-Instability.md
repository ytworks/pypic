[English](Two-Stream-Instability) · **日本語** · ソース: [`examples/two_stream`](https://github.com/ytworks/pypic/tree/main/examples/two_stream)

# 二流体不安定性とブネマン不安定性

互いに通り抜ける 2 本の粒子ビームは不安定です。ドリフト方向の小さな静電波がビームを疎密にし、その
電荷の偏りが波を強め、波はビーム粒子を捕捉するまで指数関数的に成長します。pypic はこれを 3 次元で
5 通りの設定について計算し、線形理論が厳密に予言する最大成長モードの成長率を測定します。5 つの測定値は
すべて ±10 % のゲート内にあり、ずれは {{ r.buneman.relative_error | +.1% }}（ブネマン）から
{{ r.ei_implicit.relative_error | +.1% }}（陰解法）の範囲です。

## 物理

プランナーは、全ビーム s（密度 n_s、ドリフト速度 v_s、ローレンツ因子 γ_s、熱速度 v_ts）について相対論的な
暖かい流体の分散関係を解き、最大成長波数 k_max とその成長率 γ_max を求めます。

```math
1 = \sum_s \frac{\omega_{ps}^2 / \gamma_s^3}{(\omega - k v_s)^2 - 3 k^2 v_{ts}^2},
\qquad \omega = \omega_r + i\gamma .
```

計算領域は k_max の波長 4 つ分の長さなので、最大成長モードは E_x の m = 4 フーリエモードです。線形段階では
その振幅が e^{γ_max t} で成長するはずです。

| 設定 | 構成 | 理論 k_max [ω_pe/c] | 理論 γ_max [ω_pe] | 理論 ω_r [ω_pe] |
|---|---|---|---|---|
| `ei` | 静止イオン（m_i/m_e = {{ c.ei.params.mass_ratio | .0f }}）中を u = ±{{ c.ei.params.u_drift }} c で進む 2 本の電子ビーム | {{ c.ei.planner.k_max | .3f }} | {{ c.ei.planner.gamma_max | .4f }} | 0 |
| `ep` | ペアプラズマ: 電子と陽電子がそれぞれ u = ±{{ c.ep.params.u_drift }} c の 2 本のビーム | {{ c.ep.planner.k_max | .3f }} | {{ c.ep.planner.gamma_max | .4f }} | 0 |
| `eip` | イオン（m_i/m_e = {{ c.eip.params.mass_ratio }}）上の ±{{ c.eip.params.u_drift }} c の電子ビーム、陽電子比 {{ c.eip.params.positron_fraction }} | {{ c.eip.planner.k_max | .3f }} | {{ c.eip.planner.gamma_max | .4f }} | {{ c.eip.planner.omega_r | .3f }} |
| `ei_implicit` | `ei` をエネルギー保存型の陰解法で計算 | {{ c.ei_implicit.planner.k_max | .3f }} | {{ c.ei_implicit.planner.gamma_max | .4f }} | 0 |
| `buneman` | 静止イオン中を u = {{ c.buneman.params.u_drift }} c で進む 1 本の電子ビーム | {{ c.buneman.planner.k_max | .3f }} | {{ c.buneman.planner.gamma_max | .4f }} | {{ c.buneman.planner.omega_r | .4f }} |

- **対称ビーム**（`ei`、`ep`）: 不安定根は純成長（ω_r = 0）なので、波のパターンは止まったまま成長します。
- **ブネマン不安定性**: 電子が静止したイオンの中を流れます。イオンの慣性のため成長は遅く、冷たいプラズマの
  極限では水素で γ = (√3/2)(m_e/2m_i)^{1/3} ω_pe = 0.0561 ω_pe です。プランナーは上表の暖かい根を使います。
  実周波数をもつので波は伝播します。
- **陰解法**（`ei_implicit`）: Crank–Nicolson の場とインプリシット中点法のプッシュにより、全エネルギーが
  厳密に保存されます。Δx ≈ 13 λ_De、cΔt/Δx = 1 と、陽解法の安定条件を超えた設定で計算します。

## 実験条件

規格化単位: c = ε₀ = 1、長さは c/ω_pe、時間は 1/ω_pe、密度は n₀、運動量は u = γv/c、温度は θ = T/(m c²)。

{{ table:conditions }}

<details><summary>プランナーのパラメータ（各設定の物理入力）</summary>

{{ table:params }}

</details>

## アニメーション

### 電子–イオン二流体不安定性（`ei`）

![ei 二流体不安定性のアニメーション](images/two_stream/ei.gif)

- **左上、電子の位相空間 (x, u_x)**: u_x = ±{{ c.ei.params.u_drift }} の 2 本のビームが k_max の波長で
  波打ち始め、渦（波に捕捉された粒子）に巻き上がり、最後は幅の広い台地状の分布に混ざり合います。
- **右上、E_x(x, t)**（y・z 方向に平均）: 熱雑音から 4 波長ぶんの波が成長します。ω_r = 0 なので
  パターンは動きません。
- **左下、x–y 断面の E_x**: 波面は平面です。k が x 方向の静電波だからです。
- **右下、エネルギー**: 電場エネルギーは対数軸上で直線的に（2γ で）増え、ビームが捕捉されると飽和します。
  磁場エネルギーも少し増えます。3 次元のビームは横方向のフィラメンテーションにも不安定で、プランナーはその
  成長率を γ ≈ {{ c.ei.planner.filamentation_gamma_est | .2f }} ω_pe と見積もっています。

### ブネマン不安定性（`buneman`）

![ブネマン不安定性のアニメーション](images/two_stream/buneman.gif)

成長は `ei` より約 6 倍遅くなります。E_x(x, t) の縞が傾いているのは、波が進む（ω_r ≠ 0）ためです。飽和すると
波が電子ビームを大きな位相空間の渦に捕捉し、イオンの運動エネルギーは 1 桁増えます。この不安定性は電子の
ドリフトからエネルギーを取り出し、その一部をイオンに渡します。

## 線形理論との比較

### 最大成長モードの成長

![全設定の k_max モードの成長](images/two_stream/growth.png)

各設定の k_max フーリエモードの |E_x|。破線は線形理論の成長率の傾き、灰色の帯は指数フィットの区間です。
測定された指数成長の段階は、飽和までの 2〜3 桁にわたって理論の傾きに沿っています。

### 成長率の測定値と予言値

![成長率の測定値と予言値の比](images/two_stream/growth_ratio.png)

| 設定 | 理論 γ [ω_pe] | PIC γ [ω_pe] | 誤差 | ゲート | max \|∇·E − ρ\| | max \|∇·B\| |
|---|---|---|---|---|---|---|
| `ei` | {{ r.ei.gamma_theory | .4f }} | {{ r.ei.gamma_fit_kmax | .4f }} | {{ r.ei.relative_error | +.1% }} | ±10 % | {{ r.ei.gauss_max | .1e }} | {{ r.ei.divB_max | .1e }} |
| `ep` | {{ r.ep.gamma_theory | .4f }} | {{ r.ep.gamma_fit_kmax | .4f }} | {{ r.ep.relative_error | +.1% }} | ±10 % | {{ r.ep.gauss_max | .1e }} | {{ r.ep.divB_max | .1e }} |
| `eip` | {{ r.eip.gamma_theory | .4f }} | {{ r.eip.gamma_fit_kmax | .4f }} | {{ r.eip.relative_error | +.1% }} | ±10 % | {{ r.eip.gauss_max | .1e }} | {{ r.eip.divB_max | .1e }} |
| `ei_implicit` | {{ r.ei_implicit.gamma_theory | .4f }} | {{ r.ei_implicit.gamma_fit_kmax | .4f }} | {{ r.ei_implicit.relative_error | +.1% }} | ±10 % | {{ r.ei_implicit.gauss_max | .1e }} | {{ r.ei_implicit.divB_max | .1e }} |
| `buneman` | {{ r.buneman.gamma_theory | .4f }} | {{ r.buneman.gamma_fit_kmax | .4f }} | {{ r.buneman.relative_error | +.1% }} | ±10 % | {{ r.buneman.gauss_max | .1e }} | {{ r.buneman.divB_max | .1e }} |

ガウスの法則と ∇·B は全 run で丸め誤差の範囲で成り立っています（基準はそれぞれ 1e-10 と 1e-12）。電荷保存は
電流の割り付け自体から得られ、ポアソン方程式による補正は一切使っていません。

### エネルギー保存: 陽解法と陰解法

![全エネルギー誤差、陽解法と陰解法](images/two_stream/energy.png)

陰解法は全エネルギー（粒子＋場）を run 全体で初期値の {{ r.ei_implicit.energy_error_max | .1e }} 以内に
保っており、基準の 1e-8 を大きく下回ります。陽解法では {{ r.ei.energy_error_max | .1e }} だけずれます。
陽解法の PIC はエネルギーを厳密には保存しないので、陽解法ではこの量は監視のみです。

<details><summary>各設定の診断図</summary>

各図は、エネルギー履歴、保存則の履歴、フィット付きの全 E_x モード、計算途中と終了時の電子位相空間、
E_x(x, t) を示します。

![ei](images/two_stream/ei-summary.png)
![ep](images/two_stream/ep-summary.png)
![eip](images/two_stream/eip-summary.png)
![ei_implicit](images/two_stream/ei_implicit-summary.png)
![buneman](images/two_stream/buneman-summary.png)

</details>

## 補足

- フィット区間は自動で選びます（`pypic.analysis.spectral.fit_growth`）。`eip` では陽電子が対称性を崩して
  ω_r ≠ 0 になるため、2 波モデルでフィットします。
- `buneman` では `remove_mean_current` を有効にしています。周期境界では、ドリフトする電子の正味の電流が
  一様な E_x を駆動してビーム全体を減速させてしまうためです。平均電流を差し引くことは戻り電流をモデル化
  したものに当たります。

## 再現方法

```bash
python -m pypic.launch -n 4 examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml \
    --out output/wiki/two_stream/ei
python examples/two_stream/analyze.py output/wiki/two_stream/ei --movie
```

この wiki のすべての run・図・ページは `python scripts/build_wiki.py run`、続けて `collect` と `render` で
作り直せます（[`scripts/build_wiki.py`](https://github.com/ytworks/pypic/blob/main/scripts/build_wiki.py) 参照）。
`run` が作るのはゲート規模の新規 run だけです。`collect` には他のページの記録済み本番 run の出力も必要で、
これは git に入っていません（[ホーム](Home-ja#ページの作り方) 参照）。`render` だけならチェックイン済みの
`docs/wiki/data.json` から実行できます。

## データの出所

{{ table:provenance }}

## 参考文献

- O. Buneman, "Dissipation of currents in ionized media", Phys. Rev. 115, 503 (1959).
- C. K. Birdsall and A. B. Langdon, *Plasma Physics via Computer Simulation* (IOP, 1991).
