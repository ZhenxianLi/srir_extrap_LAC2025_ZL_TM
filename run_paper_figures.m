%% RUN_PAPER_FIGURES
%  Reproduce every figure from the LAC-25 paper and write audio samples.
%  一键复现 LAC-25 论文里的全部图片，并导出参考音频。
%
%  Open this file in MATLAB and press F5 — that is the entire "one click".
%  打开本文件后按 F5（或 Run）即可，这就是论文里说的"一键复现"。
%
%  Paper / 论文:
%    Zhenxian Li and Thomas McKenzie,
%    "Spatial Room Impulse Response Extrapolation using the Image Source Method,"
%    Proc. 19th Linux Audio Conference (LAC-25), Villeurbanne, France, 2025.
%
%  Outputs / 输出:
%    figures/   — every PDF used in the paper / 论文用到的所有 PDF
%    audio_out/ — binaural-rendered reference audio / 双耳渲染的参考音频
%
%  Author / 作者: Zhenxian Li
% =========================================================================

clearvars;  close all;  clc;

% --- 0) make utils/ visible so installDependencies / setupPaths resolve ----
% --- 0) 先把 utils/ 加进 path 以便后面两个函数可见 -------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'));

% --- 1) check & auto-fetch every external resource in one call -------------
% --- 1) 一次性检查并自动获取所有外部资源 -----------------------------------
% Submodules and the 565 MB SOFA dataset are required and will be reported
% / fetched here. AMT (mckenzie2025) is optional — if you want PBC-2
% colouration reproduced, call installDependencies('InstallAMT', true) once
% from the MATLAB command window before running this script.
% 子模块与 565 MB SOFA 数据集为必需项，由本函数检查并自动获取。
% AMT（mckenzie2025）为可选项；如需复现 PBC-2 音色子图，可先在命令窗
% 单独运行 installDependencies('InstallAMT', true) 完成安装。
deps     = installDependencies();
sofaFile = deps.sofaFile;     %#ok<NASGU>  % used by srir_extrap_imgSrc_subspace_LAC

% --- 2) add every submodule to the MATLAB path ----------------------------
% --- 2) 把所有 submodule 加入 MATLAB path ---------------------------------
repoRoot = setupPaths();

% --- 3) initialise SOFA (newer SOFAtoolbox uses the same entry name) ------
% --- 3) 初始化 SOFA（新版 SOFAtoolbox 入口仍叫 SOFAstart）-----------------
SOFAstart;

% --- 4) load the binaural Ambisonic decoder used by McKenzie's PBC model ---
% --- 4) 加载论文里用到的双耳 Ambisonic 解码器 ------------------------------
% The bundled load_ambisonic_configuration script ends with a diagnostics
% tail that uses a relative `hrirs/*.wav` path and the AMT mckenzie2022
% model. The wrapper in utils/build_ambisonic_decoder.m handles both:
% it cd's to the submodule root and swallows benign diagnostic failures.
% 该 submodule 自带脚本最后会跑一段诊断（依赖相对路径 hrirs/*.wav 和
% AMT 的 mckenzie2022 模型）；wrapper 会修正 cwd 并容错地跳过诊断尾巴。
fprintf('Building binaural Ambisonic decoder / 构建双耳 Ambisonic 解码器...\n');
build_ambisonic_decoder(repoRoot);    % populates SH_ambisonic_binaural_decoder in base

% --- 5) Run the LAC reproduction script. It saves PDFs into the cwd. ------
% --- 5) 跑 LAC 复现脚本（PDF 会保存到当前工作目录）------------------------
% We cd into figures/ so all exportgraphics calls land there.
% 切到 figures/ 目录，让 exportgraphics 写出来的文件都落在那。
figDir = fullfile(repoRoot, 'figures');
if ~isfolder(figDir); mkdir(figDir); end
prevDir = pwd;
cleaner = onCleanup(@() cd(prevDir));    %#ok<NASGU>  fires on error only
cd(figDir);

fprintf('Running srir_extrap_imgSrc_subspace_LAC.m ...\n');
srir_extrap_imgSrc_subspace_LAC;

% Explicitly cd back on the success path too — `cleaner` only fires when
% the enclosing scope dies, but this script keeps living in `base` after
% `run` returns, so a chained `run_beamforming` would otherwise inherit
% the figures/ cwd.
% 显式 cd 回去 —— onCleanup 只在作用域销毁时触发，但脚本跑完仍留在 base，
% 所以后续 `run_beamforming` 会继承 figures/ 的 cwd 找不到自己。
cd(prevDir);

% --- 6) Cache LAC algorithm outputs for run_beamforming.m to pick up ------
% --- 6) 把 LAC 算法输出存盘，供 run_beamforming.m 直接加载 ----------------
% Avoids re-running the full paper pipeline just to get srir_new. The cache
% lives under figures/ (already in .gitignore).
% 这样 run_beamforming.m 不用重跑整套论文流程就能拿到 srir_new。
% 缓存文件放在 figures/ 下（已经在 .gitignore 里）。
cachePath = fullfile(figDir, 'lac_outputs.mat');
save(cachePath, 'srir_orig', 'srir_new', 'srir_tar', 'fs', ...
    'idx_LS_src_orig', 'idx_LS_rec_orig', 'idx_LS_src_tar', 'idx_LS_rec_tar');
fprintf('LAC outputs cached / LAC 输出已缓存: %s\n', cachePath);

% --- 7) Done. ---------------------------------------------------------------
fprintf('\n==============================================\n');
fprintf(' Paper figures written to / 论文图片已生成到:\n   %s\n', figDir);
fprintf(' Reference audio (if any) in / 参考音频在:\n   %s\n', fullfile(repoRoot, 'audio_out'));
fprintf('==============================================\n');
