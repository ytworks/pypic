# pypic

[English](README.md) | **日本語**

GPU で高速化した 3 次元・特殊相対論的・電磁 Particle-in-Cell (PIC) エンジンです。
ホットループは CUDA C カーネルで、CuPy の NVRTC で JIT コンパイルします (nvcc は不要)。
設定・解析は Python です。複数 GPU で動かす場合は、NCCL による 3 次元領域分割で 1 プロセスが 1 GPU を担当します。

## 特徴

| 項目 | 陽解法 `solver: explicit` (既定) | 陰解法 `solver: implicit` (エネルギー厳密保存) |
|---|---|---|
| 粒子推進 | 相対論 Buneman–Boris (`u = γv`)。磁場回転で \|u\| を厳密に保存 | 相対論 implicit midpoint。`v̄ = (u⁰+u¹)/(γ⁰+γ¹)` により `m(γ¹−γ⁰) = q Δt v̄·E` が厳密に成り立つ |
| 電磁場 | Yee 格子 FDTD (leapfrog)。∇·B = 0 を機械精度で保持 | Crank–Nicolson Maxwell。前処理付き CG で解く |
| 電流付与 | Esirkepov 法。連続の式、すなわちガウスの法則を機械精度で保存 | 軌道をセル境界で分割する Chen–Chacón 型。場の補間と電流付与が同じ重みの随伴なので、電荷と全エネルギーを同時に保存 |
| 超粒子形状 | 立方体 (CIC = 1 次 B-spline) | 同左 |
| 非線形解法 | — | 外側は場の Picard 反復。粒子ごとの軌道は Picard → Newton (信頼領域 + 直線探索) → 1 次元ブラケット解法の順に切り替える |

- **境界条件**: 軸ごと・面ごとに `periodic` / `open` / `reflecting` を指定できます。
  - `periodic`: 周期境界。
  - `open`: 陽解法では CPML、陰解法では整合導電吸収層で場を吸収します。境界を出た粒子は削除し、そのエネルギーと電荷を記録します。
  - `reflecting`: 完全導体壁と粒子の鏡面反射です。反射時は軌道を壁で分割して電流を付与するので、電荷保存は崩れません。
- **粒子種**: 電子、イオン (`mass_ratio = m_i/m_e` と価数 `Z` を指定可)、陽電子。
- **配合**: `e-i` / `e-p` / `e-i-p` から選び、常に電気的中性 (`n_e = Z n_i + n_p`) で生成します。全ビームを同じ位置に配置するので、初期の ρ は厳密に 0、E も 0 です。
- **自動パラメータ設定**: 線形分散関係から最大成長モードを求め、箱の大きさ・格子・時間刻み・計算時間・ppc・メモリ・GPU 分割を決めて、各制約に対する余裕をレポートします。どの値も上書きできます。
- **初期分布**: Maxwell–Jüttner 分布に相対論的ドリフトを与えます (Zenitani 2015)。乱数はカウンタベースなので、GPU 数を変えても同じ初期粒子になります。
- **診断**: 運動エネルギー (粒子種別)、電場・磁場エネルギー、全エネルギー収支 (境界からの流出と吸収層での散逸を含む)、運動量、max\|∇·E−ρ\|、max\|∇·B\|、E_x のモード (複素フーリエ係数)、(x, u_x) 位相空間。

単位系は正規化単位で、`c = ε0 = 1`、長さ `c/ω_pe`、時間 `1/ω_pe`、電荷 `e`、質量 `m_e`、密度 `n0` (全電子密度) です。
温度は `m_e c²` 単位 (`T_e`) で与えます。`T_e_eV` のように eV で与えることもできます。

## インストール

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[multigpu,test]"
```
CUDA 13 ドライバと CuPy (`cupy-cuda13x`) が必要です。マルチ GPU には `nvidia-nccl-cu13` が必要です。
Blackwell (sm_120) で動作を確認しています。

## 使い方: 二流体不安定 (two-stream) の例

```bash
# パラメータ計画だけを表示する
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml --plan-only
# 1 GPU で実行する
python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml
# 8 GPU で実行する
python -m pypic.launch -n 8 examples/two_stream/run.py --config examples/two_stream/configs/ep.yaml
# 解析する (線形理論との成長率比較、保存則、位相空間 → output/<name>/summary.png)
python examples/two_stream/analyze.py output/ep
```

用意している設定ファイル:

| ファイル | 配合 | 内容 |
|---|---|---|
| `configs/ei.yaml` | 電子・イオン | 対向する電子ビーム 2 本 + 静止イオン (`mass_ratio` は可変) |
| `configs/ep.yaml` | 電子・陽電子 | 中性のペアビームが対向 (`mode: symmetric`) |
| `configs/eip.yaml` | 電子・イオン・陽電子 | 正電荷の 30% を陽電子が担う (`positron_fraction`) |
| `configs/ei_implicit.yaml` | 電子・イオン | e-i と同じ物理を陰解法 (エネルギー厳密保存) で計算 |

`--set key=value` で YAML の値を上書きできます。例: `--set u_drift=0.3 --set solver=implicit`。
`SimConfig` の任意の項目は `--set overrides.ppc=64` のように上書きします。

### 主なパラメータ (`TwoStreamParams`)

- **プラズマ**
  - `mixture`: `e-i` / `e-p` / `e-i-p`
  - `mass_ratio`, `Z`, `positron_fraction`
  - `mode`
    - `electron`: 電子が対向ビーム、正イオン類は静止
    - `symmetric`: 全粒子種を中性の対向ビームにする
    - `buneman`: 電子が正イオン類に対してドリフト
  - `u_drift` (γv/c), `beam_fraction`
  - `T_e`, `T_i`, `T_p` (または `*_eV`)
  - `B0`
- **数値**
  - `solver`, `n_wavelengths`, `cells_per_wavelength`, `transverse_ratio`
  - `ppc` (位置/セル。各ビームが 1 位置に 1 粒子)
  - `cfl`, `wpe_dt_max`, `growth_times` (計算時間を `1/γ_max` 単位で)
  - `boundaries`, `max_cells`

### planner が決めるもの

1. 温かい流体の相対論的分散関係 `1 = Σ ω_ps²/γ_s³ / ((ω − k v_s)² − 3k²v_ts²)` を多成分で解き、`k_max`, `γ_max`, `ω_r` を求めます。
2. `L_x` を最大成長波長の `n_wavelengths` 倍にとり、セル数は 8 の倍数に丸めます。
3. `Δt = min(` CFL (陽解法: `cΔt/Δx = 0.95/√3`、陰解法: 1)`,` `ω_pe Δt ≤ 0.2` (陰解法 0.5)`,` `ω_ce Δt ≤ 0.2``,` `γ_max Δt ≤ 0.05``,` `|ω_r| Δt ≤ 0.2 )` とします。
4. 陽解法で `Δx > π λ_De` のときは、有限格子加熱の警告を出し、陰解法を勧めます。
5. 粒子数、GPU あたりのメモリ、rank の配置を見積もります。

## ライブラリとして使う

```python
from pypic import SimConfig, SpeciesConfig, Beam, Simulation
cfg = SimConfig(n_cells=(64, 32, 32), dx=(0.1,) * 3, dt=0.05, n_steps=1000, solver="explicit",
                boundaries=(("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                species=[SpeciesConfig("electron", -1, 1, [Beam(1.0, (0.1, 0, 0), 1e-3)]),
                         SpeciesConfig("ion", 1, 1836, [Beam(1.0, (0, 0, 0), 1e-6)])], ppc=16)
sim = Simulation(cfg)      # 複数 GPU で動かすときは python -m pypic.launch -n N で起動する
sim.run(progress=100)
sim.save()                 # output/history.npz
```

## テスト

```bash
pytest tests        # 2 GPU 以上ある場合はマルチ GPU の一致テストも実行される
```
テストの内容:
- 相対論的ジャイロ運動での \|u\| の保存と回転角、E×B ドリフト
- ガウスの法則と ∇·B (陽解法、3 種の境界)
- 陰解法の全エネルギー (相対誤差 < 1e-11、3 種の境界)
- 真空中の電磁波の Yee 数値分散、CPML の吸収 (残留エネルギー < 1e-4)
- 反射境界と開放境界での粒子の扱い
- 配合の中性
- 分散関係ソルバーと解析解の比較
- 3 次元 two-stream の成長率と理論の比較 (15% 以内)
- 1 GPU と 2 GPU の結果の一致 (相対差 < 1e-10)

## 検証結果 (8× RTX PRO 6000 Blackwell)

3 次元 two-stream の例を 8 GPU で、各設定の既定パラメータのまま実行した結果です。
γ は、k_max に対応するモードの線形成長期をフィットした値です。e-i-p は ω_r ≠ 0 なので、対向する 2 つの成長波をまとめてフィットしています。

| 設定 | 格子 / 粒子数 | γ_fit | γ_theory | 誤差 | 最大 \|ΔE_tot\|/E₀ | max\|∇·E−ρ\| | max\|∇·B\| |
|---|---|---|---|---|---|---|---|
| e-i (陽解法) | 128×32×32 / 6.3M | 0.351 | 0.342 | +2.6% | 2.2e-3 | 2.1e-12 | 4.6e-15 |
| e-p (陽解法) | 128×32×32 / 8.4M | 0.453 | 0.468 | −3.2% | 2.4e-3 | 9.2e-13 | 9.4e-15 |
| e-i-p (陽解法) | 128×32×32 / 8.4M | 0.320 | 0.332 | −3.7% | 3.5e-4 | 3.4e-12 | 8.7e-15 |
| e-i (陰解法) | 64×16×16 / 0.8M | 0.349 | 0.342 | +2.1% | **1.9e-10** | 3.4e-13 | 1.6e-15 |

- **マルチ GPU の一致**: 1 GPU と 2 GPU (2×1×1, 1×2×1, 1×1×2)、8 GPU (4×1×2) の結果は、ソルバー 2 種 × 境界 3 種の 6 通りすべてで相対差 約 1e-15 です。
- **陰解法の全エネルギー**: 周期境界、反射境界、開放境界 (吸収層での散逸を含む収支) のいずれでも 1e-13 程度です。
- **陽解法の strong scaling** (256×64×64 セル、5000 万粒子):

  | GPU 数 | 1 | 2 | 4 | 8 |
  |---|---|---|---|---|
  | 実行時間 | 11.6 s | 7.0 s | 4.1 s | 2.9 s |
  | 速度向上 | 1.0× | 1.7× | 2.8× | 4.0× |

  1 GPU では 1.0 ns/粒子/ステップです。
- **陰解法の速度**: 1 ステップは陽解法のおよそ 20 倍かかります (Picard 5〜6 回 × CG 約 15 回)。そのかわり、Δx ≫ λ_De や cΔt/Δx ≥ 1 を使えます。

## ファイル構成

```
pypic/
  kernels/common.cuh     補間、乱数、warp reduction
  kernels/fields.cu      Yee 更新 (CPML)、curl、div
  kernels/particles.cu   Boris + Esirkepov (反射壁での分割)、ρ、ロード (Maxwell-Jüttner)
  kernels/implicit.cu    陰的粒子推進 (軌道分割付与、Newton / ブラケット解法、エネルギー閉包)
  grid.py                Yee 格子、ghost 交換、境界条件、CPML、場の診断
  particles.py           SoA 粒子、ソート、rank 間移動
  simulation.py          陽解法・陰解法のメインループ、CG、診断、保存
  parallel/              NCCL 通信、3 次元領域分割
  launch.py              マルチ GPU ランチャ
  setup/                 配合 (species)、分散関係 (dispersion)、パラメータ自動設定 (planner)
examples/two_stream/     run.py, analyze.py, configs/*.yaml
tests/
```

## 注意点と制限

- **陽解法のエネルギー**: Boris 回転は \|u\| を厳密に保存しますが、系全体のエネルギーは厳密には保存しません (通常 1e-3 程度)。全エネルギーを機械精度で保存したい場合は陰解法を使ってください。陰解法の誤差は Picard の許容誤差 (`implicit.picard_tol`) 程度、およそ 1e-10 になります。
- **陰解法の非収束粒子**: ノード面や壁をかすめる粒子で個別の反復が収束しない場合でも、付与した経路に沿った仕事に一致するよう \|u¹\| を補正します。このため全エネルギーと電荷の保存は保たれます。非収束の件数は `sim.nonconverged` で確認できます。
- **陽解法の開放境界**: エネルギー収支のうち、ポインティング流束の評価は近似です。
- **陰解法の開放境界**: 整合導電層を使います。散逸エネルギーは厳密に計上されますが、吸収性能は CPML より劣ります。
- **`mode: buneman`**: 正味の電流を持つので、k=0 モードが ω_pe で振動します。
- **精度**: 現在は倍精度のみです。今後の候補は、float32 モード、共有メモリを使った電流付与、通信と計算のオーバーラップです。
