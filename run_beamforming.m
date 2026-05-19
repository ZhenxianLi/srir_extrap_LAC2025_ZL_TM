%% RUN_BEAMFORMING
%  Three-panel comparison of the LAC paper method vs. the beamforming-
%  refined extension (LAC + beam_filter inside each arrival loop), using
%  the exact same paper LAC pipeline twice — once with the beamforming
%  flag off, once with it on.
%  论文 LAC 方法 vs 波束成形细化扩展（在 LAC 主循环每个反射内部插入
%  beam_filter）的三栏对比。两次都跑同一份论文 LAC 流程，只切换
%  flag_use_beamforming 开关。
%
%  Panel layout matches Tom McKenzie's email screenshot:
%    Original SRIR  |  LAC extrapolation  |  Beamforming-refined
%  三栏布局与 Tom 邮件截图一致：原始 SRIR / LAC 外推 / 波束成形细化。
%
%  Pre-requisite: run `run_paper_figures.m` first so the LAC outputs are
%  cached in figures/lac_outputs.mat. This script then loads the cache,
%  re-runs only the LAC algorithm with flag_use_beamforming=true (silent
%  — no figures, no audio), and renders the three-panel comparison.
%  前置：先跑一次 run_paper_figures.m 生成 figures/lac_outputs.mat 缓存。
%  本脚本加载缓存，再以静默方式（不画图、不写音频）跑一次 LAC 算法
%  （flag_use_beamforming=true），最后画三栏对比。
%
%  Author / 作者: Zhenxian Li
% =========================================================================

clearvars; close all; clc;

% --- 0) path ----------------------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'));

% --- 1) dependencies + paths + SOFA -----------------------------------------
installDependencies();
repoRoot = setupPaths();
SOFAstart;

% --- 2) Load cached LAC outputs from run_paper_figures.m --------------------
cachePath = fullfile(repoRoot, 'figures', 'lac_outputs.mat');
if ~exist(cachePath, 'file')
    error('run_beamforming:cacheMissing', ...
        ['LAC output cache not found at:\n  %s\n' ...
         'Please run  run_paper_figures.m  first.\n' ...
         '未找到 LAC 输出缓存：\n  %s\n请先运行 run_paper_figures.m。'], ...
         cachePath, cachePath);
end
fprintf('Loading cached LAC outputs / 加载 LAC 输出缓存:\n  %s\n', cachePath);
cache = load(cachePath);
srir_orig    = cache.srir_orig;
srir_lac     = cache.srir_new;    % LAC paper method / 论文 LAC 方法
srir_tar     = cache.srir_tar;    % ground-truth target / 目标真值
fs           = cache.fs;

% --- 3) Re-run the paper LAC script with beamforming flag ON, silently ----
% --- 3) 用 flag_use_beamforming=true 静默重跑论文 LAC 脚本 -----------------
fprintf('\nRunning LAC pipeline again with flag_use_beamforming = true ...\n');

% Build the decoder again (paper script needs SH_ambisonic_binaural_decoder
% in base for binaural-rendering steps even though we suppress audio).
% 重新构建 decoder（论文脚本需要它，即便我们关掉了音频导出）。
build_ambisonic_decoder(repoRoot);

% Push the LAC inputs + beamforming flag + silence flags into base so the
% script picks them up. The script does NOT clear locals (we removed
% `clear all` in our reproducible cleanup), so it will see these vars.
% 把 LAC 输入参数 + beamforming flag + 静默 flag 推到 base 工作空间，
% 让脚本能识别。脚本里不会 clear，能读到这些变量。
assignin('base', 'idx_LS_src_orig',      cache.idx_LS_src_orig);
assignin('base', 'idx_LS_rec_orig',      cache.idx_LS_rec_orig);
assignin('base', 'idx_LS_src_tar',       cache.idx_LS_src_tar);
assignin('base', 'idx_LS_rec_tar',       cache.idx_LS_rec_tar);
assignin('base', 'flag_use_beamforming', true);
assignin('base', 'beam_beta',            0.8);
assignin('base', 'beam_sigma',           1.2);
assignin('base', 'save_figs',            false);   % no PDF writes
assignin('base', 'flag_render_binaural', false);   % no audio writes
assignin('base', 'repoRoot',             repoRoot);

% Hide figures the paper script keeps creating along the way.
% 关掉论文脚本一路创建的所有 figure 窗口（batch 下本来就不显示，加这行
% 让交互模式也安静）。
prevVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
visRestore = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));   %#ok<NASGU>

% Run the paper LAC script as a script (so its variables land in base).
% It will produce a fresh srir_new — this time the beamforming-refined one.
% 把论文脚本当作脚本跑（变量落到 base），生成新的 srir_new —— 这次是
% 波束成形细化版本。
evalin('base', 'srir_extrap_imgSrc_subspace_LAC');
srir_beam = evalin('base', 'srir_new');

% Close all the (hidden) figures the paper script created.
% 关掉脚本创建的（隐藏）figure 窗口。
close all;
clear visRestore;        % triggers visibility restore

fprintf('Beamforming-refined SRIR obtained / 波束成形细化结果已得到\n');

% --- 4) Three-panel spherical power comparison -----------------------------
% --- 4) 三栏球面功率对比 ---------------------------------------------------
% Use exactly the same colormap + clim as the paper's DoA plots (line 398
% of srir_extrap_imgSrc_subspace_LAC.m), so the visual style matches.
% 与论文 DoA 球面图一致：flipud(bone) + clim([0 3])。
figDir = fullfile(repoRoot, 'figures');
if ~isfolder(figDir); mkdir(figDir); end

figh = figure('Color', 'w', 'Position', [60 60 1900 380], 'Visible', 'on');

subplot(1, 4, 1);
[P, ~, ~, gd] = get_pwd(srir_orig, fs);
heatmap_plot(rad2deg(gd(:,1)), rad2deg(gd(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
title('Original SRIR');

subplot(1, 4, 2);
[P, ~, ~, gd] = get_pwd(srir_lac, fs);
heatmap_plot(rad2deg(gd(:,1)), rad2deg(gd(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
title('LAC extrapolation');

subplot(1, 4, 3);
[P, ~, ~, gd] = get_pwd(srir_beam, fs);
heatmap_plot(rad2deg(gd(:,1)), rad2deg(gd(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
title(sprintf('Beamforming-refined (\\beta=%.2g, \\sigma=%.2g)', ...
    evalin('base', 'beam_beta'), evalin('base', 'beam_sigma')));

subplot(1, 4, 4);
[P, ~, ~, gd] = get_pwd(srir_tar, fs);
heatmap_plot(rad2deg(gd(:,1)), rad2deg(gd(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
title('Target (ground truth)');

cb = colorbar('Position', [0.945 0.18 0.010 0.7]);
ylabel(cb, 'Normalised power (dB)');

outPdf = fullfile(figDir, 'beamforming_4panel.pdf');
outPng = fullfile(figDir, 'beamforming_4panel.png');
exportgraphics(figh, outPdf);
exportgraphics(figh, outPng, 'Resolution', 200);

fprintf('\n==============================================\n');
fprintf(' Four-panel comparison saved / 四栏对比图已保存:\n');
fprintf('   %s\n', outPdf);
fprintf('   %s\n', outPng);
fprintf('==============================================\n');

% =========================================================================
function [P_pwd, doa_est, doa_est_P, grid_dirs] = get_pwd(srir, fs)
% Local copy of get_pwd from srir_extrap_imgSrc_subspace_LAC.m (line ~653)
% so this script can plot using exactly the same SRP+peak-pick definition
% as the paper's Fig. 6.
% 论文脚本里 line 653 那个 get_pwd 的本地副本：让 run_beamforming.m 用与
% 论文图 6 完全相同的 SRP+峰值定义。
degreeResolution    = 2;
order               = sqrt(size(srir, 2)) - 1;
nSrc                = 7;
numSamps            = fs * 0.3;
highPassFilterFreq  = 3000;
kappa               = 40;
grid_dirs           = grid2dirs(degreeResolution, degreeResolution, 0, 0);
P_src               = diag(ones(numSamps, 1));
[~, filtHi, ~]      = ambisonic_crossover(highPassFilterFreq, fs);
doa_est             = zeros(nSrc, 2);
doa_est_P           = zeros(nSrc, 1);
P_pwd               = zeros(size(grid_dirs, 1), 1);
Y_src               = filter(filtHi, 1, srir(1:numSamps, :));
stVec               = Y_src';
sphCOV              = stVec * P_src * stVec' + eye((order + 1)^2) / (4 * pi);
[P_pwd(:), est_dirs_pwd, est_dirs_P] = sphPWDmap_compat(sphCOV, grid_dirs, nSrc, kappa);
est_dirs_pwd        = est_dirs_pwd * 180 / pi;
est_dirs_pwd(est_dirs_pwd < -170) = est_dirs_pwd(est_dirs_pwd < -170) + 360;
doa_est(:, :)       = est_dirs_pwd;
doa_est_P(:)        = est_dirs_P;
end
