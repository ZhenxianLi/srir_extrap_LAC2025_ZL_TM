# SRIR Extrapolation via Image-Source Modelling

*Reference MATLAB implementation accompanying* **"Spatial Room Impulse Response
Extrapolation using the Image-Source Method"** *(LAC-25, 2025)*, plus the
post-paper **beamforming-refined** follow-up method.

> 单测量空间房间脉冲响应（SRIR）的镜像源外推法 —— LAC-25 论文的官方
> MATLAB 复现代码，外加论文之后做的"波束成形细化"工作。

---

## What's here / 仓库内容

| Path | Purpose |
|------|---------|
| [`run_paper_figures.m`](run_paper_figures.m) | **One-click** reproduction of every figure & audio sample in the LAC-25 paper.<br>**一键**复现 LAC-25 论文里的所有图片和参考音频。 |
| [`run_beamforming.m`](run_beamforming.m) | Beamforming-refined extrapolation; produces the 3-panel sphere comparison (Original / LAC / new) and the RMS / Colouration / Localisation table.<br>波束成形细化方法；产出三栏球面对比图，并打印 RMS / 音色 / 定位差异表。 |
| [`FinalCode/`](FinalCode/) | Core algorithm scripts (`srir_extrap_imgSrc_subspace_LAC.m`, `srir_extrap_beamforming.m`, `srirSubspaceDecomp.m`, `steerAmbiSOFA*.m`).<br>核心算法脚本。 |
| [`utils/`](utils/) | Project-local helpers (`installDependencies`, `setupPaths`, `build_ambisonic_decoder`, `sphPWDmap_compat`, `downloadSrirDataset`, `fig`, `freqplot_smooth`) and `utils/conventions/` (legacy SOFA convention CSVs patched back into SOFAtoolbox at install time).<br>项目自带的小工具，`utils/conventions/` 是装载时回填到 SOFAtoolbox 的老版 SOFA convention 定义。 |
| [`external/`](external/) | Git submodules → upstream dependencies (see below).<br>Git 子模块 → 上游依赖。 |
| [`audio_samples/`](audio_samples/) | Dry source signals used for binaural rendering.<br>用于双耳渲染的干声样本。 |
| [`6dof_SRIRs_eigenmike_SH/`](6dof_SRIRs_eigenmike_SH/) | Empty placeholder — the 565 MB SOFA dataset lands here on first run.<br>占位文件夹，首次运行会自动下载 565 MB 的 SOFA 数据集。 |
| `figures/`, `audio_out/` | Generated artefacts (in `.gitignore`).<br>运行后生成的图片和音频（不入库）。 |

---

## Clone (with submodules!) / 克隆（务必带 submodule）

The four open-source toolboxes that used to be bundled into this repo are
now **git submodules** living in `external/`. After cloning you must pull
them in, otherwise nothing will run.

之前直接打包在仓库里的四个第三方工具箱，现在都改为 **git submodule** 放在
`external/`。所以 clone 完一定要再拉一次 submodule，否则跑不起来。

```bash
# Either / 二选一:
git clone --recurse-submodules https://github.com/ZhenxianLi/srir_extrap_LAC2025_ZL_TM.git

# or, if you already cloned without --recurse-submodules:
# 或者，已经 clone 过但忘了 --recurse-submodules：
cd srir_extrap_LAC2025_ZL_TM
git submodule update --init --recursive
```

### Submodules pinned in `external/`

| Submodule | Upstream | What it provides |
|-----------|----------|------------------|
| `SOFAtoolbox` | https://github.com/sofacoustics/SOFAtoolbox | SOFA file I/O (`SOFAstart`, `SOFAload`, `SOFAplotGeometry`). Successor to the archived `API_MO`. |
| `Higher-Order-Ambisonics` | https://github.com/polarch/Higher-Order-Ambisonics | HOA encode/decode/rotate (`rotateHOA_N3D`, `getMaxREweights`, …). |
| `Spherical-Harmonic-Transform` | https://github.com/polarch/Spherical-Harmonic-Transform | SHT primitives (`grid2dirs`, `getRSH`, …). |
| `Spherical-Array-Processing` | https://github.com/polarch/Spherical-Array-Processing | Plane-wave decomposition (`sphPWDmap`) used in the DoA plots. |
| `binaural_ambisonic_preprocessing` | https://github.com/thomas-mckenzie/binaural_ambisonic_preprocessing | Binaural Ambisonic decoder & pre-processing pipeline. |

---

## Dependency installer / 依赖安装器

Everything external this project needs — git submodules, the 565 MB SOFA
dataset, and the optional AMT toolbox — is checked and (where possible)
auto-downloaded by a single helper:

本项目所有外部依赖（git 子模块、565 MB 的 SOFA 数据集、可选的 AMT 工具箱）
都由同一个脚本集中检查并按需自动下载：

```matlab
addpath('utils')          % once per MATLAB session / 每个会话一次

installDependencies()                           % default: auto-fetch SOFA
                                                % 默认：自动下 SOFA
installDependencies('InstallAMT', true)         % also install AMT (~14 MB zip)
                                                % 顺便装 AMT
installDependencies('Interactive', false)       % silent mode (CI)
                                                % 静默模式（CI 用）
```

The two main scripts (`run_paper_figures.m`, `run_beamforming.m`) call
`installDependencies()` themselves on the first line, so for the common case
**you do not need to call it manually** — just press F5 on either main
script. Run it manually only if you want to pre-install AMT or do a
silent CI-style install.

两个主脚本（`run_paper_figures.m`、`run_beamforming.m`）会在内部自动调用
`installDependencies()`，所以日常使用**不需要手动运行**，按 F5 直接跑主脚本
即可。只有当你想预装 AMT、或在 CI 上做静默安装时才需要手动调用。

---

## Quick start / 快速上手

### 1. Reproduce the paper / 复现论文图

Open MATLAB in the repo root and press **F5** on `run_paper_figures.m`.

在仓库根目录打开 MATLAB，运行 [`run_paper_figures.m`](run_paper_figures.m) 即可：

```matlab
run('run_paper_figures.m')
```

What this does step-by-step / 这一步会做这些事:

1. Call `installDependencies()` — checks every submodule, auto-downloads
   the 565 MB Arni 6DoF SOFA dataset
   (mirror: <https://zenodo.org/records/5720724>) if missing.
   调用 `installDependencies()` —— 检查所有子模块；首次运行会自动下载
   565 MB 的 Arni 6DoF SOFA 数据集。
2. Add every dependency to the MATLAB path via `setupPaths.m`.
   通过 `setupPaths.m` 把所有依赖加进 path。
3. Build the binaural Ambisonic decoder.
   构建双耳 Ambisonic 解码器。
4. Run `FinalCode/srir_extrap_imgSrc_subspace_LAC.m`, which writes every
   paper figure into `figures/` and reference audio into `audio_out/`.
   运行 `FinalCode/srir_extrap_imgSrc_subspace_LAC.m`，把论文里所有图
   存到 `figures/`，参考音频存到 `audio_out/`。

Generated artefacts / 生成的产物:

| File | Paper figure |
|------|--------------|
| `figures/srir_extrap_imgsrc_orig.pdf` | Fig. 2(a) — original 2nd-order image-source plot |
| `figures/srir_extrap_imgsrc_tar.pdf`  | Fig. 2(b) — target 2nd-order image-source plot |
| `figures/srir_extrap_TD_comparison_orig.pdf` | Fig. 3(a) — time-domain Original vs. Target |
| `figures/srir_extrap_TD_comparison_extrap.pdf` | Fig. 3(b) — time-domain Extrapolated vs. Target |
| `figures/srir_extrap_FD_comparison_orig_bin.pdf` | Fig. 4(a) — binaural frequency response, Original |
| `figures/srir_extrap_FD_comparison_tar_bin.pdf` | Fig. 4(b) — binaural frequency response, Extrapolated |
| `figures/srir_extrap_hor_doa.pdf` | Fig. 5 — horizontal DoA |
| `figures/srir_extrap_doa_{original,extrap,target}.pdf` | Fig. 6 — spherical DoA heatmaps |

### 2. Run the steerAmbiSOFA demo / 跑 steerAmbiSOFA 演示

Prerequisite: `run_paper_figures.m` must have been run at least once
(it caches `srir_orig / srir_new / srir_tar` for the demo to pick up).

前置条件：至少跑过一次 `run_paper_figures.m`（它会缓存
`srir_orig / srir_new / srir_tar`，本 demo 直接加载）。

```matlab
run('run_beamforming.m')
```

> **Status note** — this script does **NOT** reproduce the three-panel
> screenshot from Tom McKenzie's email (Original / LAC / Beamforming-
> refined). That figure requires Tom's full beamforming-refined pipeline
> ("I got the system working now quite nicely!"), which lives in his
> local workspace and has not been shared back. The script here is a
> **component demo** of [`FinalCode/steerAmbiSOFA.m`](FinalCode/steerAmbiSOFA.m)
> only — the HOA beam-steering function that was sent to Tom for
> integration. [`FinalCode/steerAmbiSOFA_maxre.m`](FinalCode/steerAmbiSOFA_maxre.m)
> exists as an exploration of max-rE weighting but did **not** make it
> into the final workflow.
>
> **状态说明** —— 本脚本**不复现** Tom 邮件截图里的三栏对比（原始 / LAC /
> Beamforming-refined）。那张图需要 Tom 完整的 beamforming-refined pipeline
> （"I got the system working now quite nicely!"），目前在他本机，尚未发回。
> 本脚本只是 [`FinalCode/steerAmbiSOFA.m`](FinalCode/steerAmbiSOFA.m) 的
> **组件演示** —— 当时发给 Tom 的 HOA 波束 steering 函数。
> [`FinalCode/steerAmbiSOFA_maxre.m`](FinalCode/steerAmbiSOFA_maxre.m) 是
> 探索 max-rE 加权的尝试，**未**进入最终工作流程。

What this demo does / 本 demo 实际做的事：

1. Load the cached `srir_orig / srir_new / srir_tar` from
   `figures/lac_outputs.mat`.
   加载 `figures/lac_outputs.mat` 里的 LAC 输出。
2. Wrap `srir_new` (LAC-extrapolated) into a 1-measurement temporary
   SOFA file and pass it through the real `steerAmbiSOFA(...)` (file in,
   file out — exercising the actual component end-to-end).
   把 `srir_new`（LAC 外推结果）包装成 1 条 IR 的临时 SOFA 文件，喂给
   真实的 `steerAmbiSOFA(...)`（文件输入 → 文件输出，走完整流程）。
3. Plot a three-panel spherical-power comparison:
   *Original / LAC-extrapolated / LAC + steerAmbiSOFA*.
   画三栏球面功率图：原始 / LAC 外推 / LAC + steerAmbiSOFA。

Outputs / 产物:

* `figures/steerAmbiSOFA_demo_3panel.pdf` / `.png` — the three-panel
  comparison described above.
  上述三栏对比图。

---

## Compatibility shims / 兼容性补丁

Several upstream toolboxes have moved on while this paper code stayed on
2021-era APIs. `installDependencies` patches the gaps for you — these
are listed here just so you know what is happening behind the scenes:

几个上游工具箱已经更新，而本论文代码停留在 2021 年的 API。
`installDependencies` 会自动 patch 这些差异，下面列出来只是让你知道
背后发生了什么：

* **SOFA convention `SingleRoomDRIR`** — the Arni 6DoF dataset declares
  this convention, but SOFAtoolbox 2.x removed it (renamed to
  `SingleRoomSRIR`). [`utils/conventions/SingleRoomDRIR_0.3.csv`](utils/conventions/SingleRoomDRIR_0.3.csv)
  is copied into the SOFAtoolbox submodule on first run; `SOFAstart`
  compiles it to `.mat`. The submodule's local `.git/info/exclude` is
  also updated so the patched CSV does not pollute `git status`.<br>
  Arni 6DoF 数据集仍使用老 convention `SingleRoomDRIR`，但 SOFAtoolbox 2.x
  已把它改名为 `SingleRoomSRIR` 并移除。首次运行时把 CSV 拷到 submodule
  里，由 `SOFAstart` 自动编译为 `.mat`，并把 submodule 的本地
  `.git/info/exclude` 也更新好，避免污染 `git status`。
* **`sphPWDmap`** (Spherical-Array-Processing) — upstream simplified to
  3 inputs / 2 outputs, dropping the user-supplied `kappa` and the
  per-peak power vector. [`utils/sphPWDmap_compat.m`](utils/sphPWDmap_compat.m)
  restores the old 4-in / 3-out behaviour by re-doing peak picking with
  the caller's `kappa` and looking up per-peak powers.<br>
  上游 `sphPWDmap` 被简化成 3 入 2 出，丢掉了 `kappa` 参数和每峰功率
  返回值。`utils/sphPWDmap_compat.m` 自己重做 peak picking 还原老 API。
* **`load_ambisonic_configuration` diagnostics tail** — the bundled
  script ends with a call to `test_ambisonic_decoder` that uses a
  relative `hrirs/*.wav` path and AMT's `mckenzie2022`. [`utils/build_ambisonic_decoder.m`](utils/build_ambisonic_decoder.m)
  cd's into the submodule root, snapshots the base workspace (the
  script does `clear variables` internally), and swallows benign
  diagnostic-tail failures.<br>
  binaural_ambisonic_preprocessing 自带的 `test_ambisonic_decoder` 用
  相对路径取 hrirs 并调用 AMT 的 `mckenzie2022`。wrapper
  `utils/build_ambisonic_decoder.m` 帮我们处理 cwd、保护 base 工作空间
  （脚本里有 `clear variables`），并优雅地接住诊断尾巴的失败。

---

## Optional: AMT toolbox & PBC-2 / 选装：AMT 工具箱与 PBC-2

The paper reports binaural colouration with the **PBC-2** model
(`mckenzie2025`) from the [Auditory Modeling Toolbox (AMT)](https://amtoolbox.org).

**Status note** — as of AMT **1.6.0** (latest release, Oct 2024),
`mckenzie2025` is **NOT yet shipped upstream**; only `mckenzie2021` and
`mckenzie2022` are. Installing AMT therefore enables every part of the
project *except* the PBC-2 plot, which still needs an `mckenzie2025.m`
obtained directly from Thomas McKenzie. Drop that file into
`external/amtoolbox/<unzipped-dir>/models/` and re-run.

**当前状态** —— 截至 AMT **1.6.0**（最新版，2024-10），`mckenzie2025`
**尚未合入主线**，仓库里只有 `mckenzie2021` 和 `mckenzie2022`。所以装
AMT 之后项目其它部分都能跑，**唯独 PBC-2 子图**还需要从 Thomas McKenzie
处单独获取 `mckenzie2025.m`，放入 `external/amtoolbox/<解压目录>/models/`
后重新运行。

Installing AMT (still recommended — `mckenzie2022` is used by the
`binaural_ambisonic_preprocessing` diagnostics):

安装 AMT 的方法（仍然推荐安装 —— 子模块 `binaural_ambisonic_preprocessing`
的诊断脚本会用到 `mckenzie2022`）：

```matlab
addpath('utils')
installDependencies('InstallAMT', true)                          % default location
installDependencies('InstallAMT', true, 'AMTRoot', ...           % custom location
    fullfile(pwd, 'external', 'amtoolbox'))
```

Default install location: `C:\MATLAB\amtoolbox` on Windows or
`~/MATLAB/amtoolbox` on Linux/macOS. The downloader fetches the small
(~14 MB) source zip from SourceForge; `amt_start` then pulls auxdata on
demand.

默认安装位置：Windows 是 `C:\MATLAB\amtoolbox`；Linux/macOS 是
`~/MATLAB/amtoolbox`。下载器会从 SourceForge 拉 ~14 MB 的源码 zip；
`amt_start` 会在需要时按需获取 auxdata。

If AMT is not available, `run_paper_figures.m` skips the PBC-2 block and
emits a warning — everything else (figures and audio) still runs.

如果没装 AMT，`run_paper_figures.m` 只会跳过 PBC-2 那一段并给出 warning，
其它流程（图片 + 音频）照常输出。

---

## MATLAB requirements / MATLAB 环境

* MATLAB R2021b or newer (the `arguments` block in `srirSubspaceDecomp.m`
  needs ≥ R2019b; `clim()` needs ≥ R2022a — use `caxis()` on older MATLAB).
  R2021b 及以上即可。
* Toolboxes used by binaural rendering:
  双耳渲染会用到这些 toolbox：
  * Signal Processing Toolbox
  * Phased Array System Toolbox (used by `binaural_ambisonic_preprocessing`)
  * Curve Fitting Toolbox (same)

---

## Citation / 引用

```bibtex
@inproceedings{li2025srir,
  author    = {Zhenxian Li and Thomas McKenzie},
  title     = {Spatial Room Impulse Response Extrapolation using the Image-Source Method},
  booktitle = {Proc. 19th Linux Audio Conference (LAC-25)},
  address   = {Villeurbanne, France},
  year      = {2025},
  month     = jun,
}
```

The 6DoF SRIR dataset:
6DoF SRIR 数据集来源：

```bibtex
@misc{mckenzie2021dataset,
  author    = {Thomas McKenzie and Leo McCormack and Christoph Hold},
  title     = {Dataset of spatial room impulse responses in a variable
               acoustics room for six degrees-of-freedom rendering and analysis},
  year      = {2021},
  doi       = {10.5281/zenodo.5720724},
  url       = {https://zenodo.org/records/5720724},
}
```

---

