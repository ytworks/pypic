[English](Particle-Orbits) · **日本語** · ソース: [`examples/particle_orbits`](https://github.com/ytworks/pypic/tree/main/examples/particle_orbits)

# 与えられた場の中のテスト粒子軌道

粒子プッシャーは PIC コードの心臓部です。ここでは、相対論的な軌道が閉じた式で分かっている場の中でテスト電子を
動かします。テスト粒子は場を感じますが電流を割り付けないので、場は与えたとおりに保たれます。すべての軌道を
毎ステップ追跡し、厳密解と比べます。どのケースでも位置誤差の最大値は軌道の大きさの
{{ r.gyration.max_rel_position_error | .2% }}（ジャイロ運動）以下で、Δt を半分にすると誤差は
{{ r.plane_wave.convergence_ratio | .2f }}（平面波）〜 {{ r.gyration.convergence_ratio | .2f }} 分の 1 になります。
スキームは 2 次精度です。

## 物理

| 設定 | 場 | 厳密解 |
|---|---|---|
| `gyration` | 一様な B_z、u⊥ = 0.1 … 10 | 角振動数 qB/(γm)、半径 u⊥m/(qB) の回転 |
| `exb` | E_y と B_z、E/B = 0.5 | E×B/B² でドリフトする案内中心のまわりのジャイロ運動。ドリフト系で厳密に解く |
| `hyperbolic` | 一様な E_x、静止から | u = qEt/m、x = (m/qE)(√(1 + (qEt/m)²) − 1) |
| `e_gt_b` | E_y > B_z（E/B > 1） | B/E で動く系での双曲線運動を実験室系に戻したもの |
| `plane_wave` | 真空中の平面波、a₀ = {{ c.plane_wave.params.a0 }}、Yee 格子の厳密な固有モード | 光面不変量 γ − u_x = 1: 8 の字軌道、ドリフト a₀²/(4 + a₀²) = {{ c.plane_wave.planner.drift_theory }} c |

解は `examples/particle_orbits/theory.py` にあり、それ自体を ODE の数値積分と 1e-11 の精度で照合しています。
相対論的ボリスプッシュは B による運動量の回転を厳密に行うので、純粋な磁場中では |u| が丸め誤差の範囲で
保存されます。

## 実験条件

規格化単位: c = 1、長さは c/ω、時間は 1/ω。ω は各ケースの特性周波数（サイクロトロン周波数または波の
周波数）です。各ケースでテスト電子 8 個、GPU 1 枚を使い、Δt と Δt/2 で 2 回ずつ計算します（平面波では Δx も
半分にします）。

{{ table:conditions }}

<details><summary>プランナーのパラメータ</summary>

{{ table:params }}

</details>

## アニメーション

![テスト粒子軌道](images/particle_orbits/orbits.gif)

各パネルでは、PIC の軌道（色付きの実線）を厳密な軌道（黒の破線）の上に描いていきます。両者は 1 % の何分の
1 という精度で一致するので、破線は PIC の線に隠れて見えません。

- **ジャイロ運動**: u⊥ = 2.7、5.2、10 の 3 個の電子。運動量が大きいほど半径が大きくなります。
- **E×B ドリフト**: x 方向にドリフトするサイクロイド。
- **双曲線運動**: 一定の E_x の下での x(t)。速度は c に近づいていきます。
- **E > B**: 電場が勝ち、電子は走り去ります。
- **平面波**: 平均ドリフトとともに動く系で見ると、軌道は古典的な 8 の字になります。

## 厳密解との比較

### 2 次の収束

![時間刻みに対する誤差](images/particle_orbits/convergence.png)

各ケースの Δt と Δt/2 での最大位置誤差（軌道の大きさで規格化）。破線は傾き 2 の目安線で、2 次精度の
スキームなら Δt を半分にすると誤差は 4 分の 1 になります。

| 設定 | Δt | Δt での最大誤差 | Δt/2 での最大誤差 | 比 | ゲート |
|---|---|---|---|---|---|
| `gyration` | {{ r.gyration.dt | .4g }} | {{ r.gyration.max_rel_position_error | .2e }} | {{ r.gyration.refined_max_error | .2e }} | {{ r.gyration.convergence_ratio | .3f }} | 誤差 ≤ 2 %、比 3〜5 |
| `exb` | {{ r.exb.dt | .4g }} | {{ r.exb.max_rel_position_error | .2e }} | {{ r.exb.refined_max_error | .2e }} | {{ r.exb.convergence_ratio | .3f }} | 誤差 ≤ 2 %、比 3〜5 |
| `hyperbolic` | {{ r.hyperbolic.dt | .4g }} | {{ r.hyperbolic.max_rel_position_error | .2e }} | {{ r.hyperbolic.refined_max_error | .2e }} | {{ r.hyperbolic.convergence_ratio | .3f }} | 誤差 ≤ 2 %、比 3〜5 |
| `e_gt_b` | {{ r.e_gt_b.dt | .4g }} | {{ r.e_gt_b.max_rel_position_error | .2e }} | {{ r.e_gt_b.refined_max_error | .2e }} | {{ r.e_gt_b.convergence_ratio | .3f }} | 誤差 ≤ 2 %、比 3〜5 |
| `plane_wave` | {{ r.plane_wave.dt | .4g }} | {{ r.plane_wave.max_rel_position_error | .2e }} | {{ r.plane_wave.refined_max_error | .2e }} | {{ r.plane_wave.convergence_ratio | .3f }} | 誤差 ≤ 2 %、比 3〜5 |

平面波のドリフト速度: 測定 {{ r.plane_wave.drift_measured | .4f }} c、理論 {{ r.plane_wave.drift_theory_node | .4f }} c。

### 各ケースの軌道と誤差

左: x–y 面の軌道。PIC（色）を厳密解（破線）に重ねたもの。右: |x_PIC − x_厳密| の時間変化。

![gyration](images/particle_orbits/gyration-orbits.png)
![E×B](images/particle_orbits/exb-orbits.png)
![hyperbolic](images/particle_orbits/hyperbolic-orbits.png)
![E > B](images/particle_orbits/e_gt_b-orbits.png)
![plane wave](images/particle_orbits/plane_wave-orbits.png)

## 補足

- **比が 2 ではなく 4 になる理由**: 陽解法は運動量を半ステップ過去に戻すところから始めるので、位置と運動量は
  最初から、リープフロッグ法の前提どおり Δt/2 ずれています。そのため全体の誤差が 2 次になります。
- **平面波**: 波は Yee 格子の固有モードなので、この run では Δt と一緒に Δx も細かくしています。その比
  （{{ r.plane_wave.convergence_ratio | .2f }}）は空間と時間をあわせた収束を表し、4 より少し小さくなります。

## 再現方法

```bash
python -m pypic.launch -n 1 examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml \
    --out output/wiki/particle_orbits/plane_wave
python -m pypic.launch -n 1 examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml \
    --set dt_scale=0.5 --out output/wiki/particle_orbits/plane_wave_half
python examples/particle_orbits/analyze.py output/wiki/particle_orbits/plane_wave \
    --refined output/wiki/particle_orbits/plane_wave_half
```

## データの出所

{{ table:provenance }}

## 参考文献

- J. P. Boris, "Relativistic plasma simulation — optimization of a hybrid code", Proc. 4th Conf. Num.
  Sim. Plasmas, 3 (1970).
- L. D. Landau and E. M. Lifshitz, *The Classical Theory of Fields*（一様な定常場と平面波の中の荷電粒子の運動）.
