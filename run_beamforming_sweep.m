%% RUN_BEAMFORMING_SWEEP
%  Coarse 5x5 grid search over (beam_beta, beam_sigma) followed by, if
%  needed, an automatic fminsearch refinement. The objective is a
%  composite loss combining three terms vs. the ground-truth target SRIR:
%
%    L = mean(L_time/L_time_LAC, L_spec/L_spec_LAC, L_pwd/L_pwd_LAC)
%
%   * L_time : W-channel time-domain NRMSE                 (energy/envelope)
%   * L_spec : 1/3-oct W-channel log-spectrum mean abs err (colouration)
%   * L_pwd  : normalised PWD spherical map L1 distance    (DoA / spatial)
%
%   Each sub-loss is divided by the LAC baseline so all three live on a
%   comparable scale, and a final composite < 1.0 means "better than LAC".
%
%   对 (beam_beta, beam_sigma) 做 5x5 粗扫；若 best 没有明显优于 LAC，
%   则用 fminsearch 在 best 附近做精细优化。损失为三项之和（时域 NRMSE、
%   1/3 oct 频谱 dB 差、球面 PWD L1 差），每项除以 LAC 基线后等权相加，
%   < 1.0 即"优于 LAC"。
%
%  Side-effects: writes
%    figures/beamforming_sweep_grid.pdf  (5x5 sphere-power grid)
%    figures/beamforming_sweep_grid.png
%    figures/beamforming_sweep_losses.csv
%    figures/beamforming_final_4panel.pdf (Original / LAC / FINAL / Target)
%    figures/beamforming_final_4panel.png
%  Does NOT touch run_paper_figures.m outputs (paper PDFs, LAC audio).
%  不影响 run_paper_figures.m 的产物（论文 PDF、LAC 音频）。
%
%  Pre-requisite: run run_paper_figures.m at least once to populate
%  figures/lac_outputs.mat.
%
%  Author / 作者: Zhenxian Li
% =========================================================================

clearvars; close all; clc;

% --- 0) path + deps + SOFA -------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'));
installDependencies();
repoRoot = setupPaths();
SOFAstart;

% --- 1) load LAC outputs cache --------------------------------------------
cachePath = fullfile(repoRoot, 'figures', 'lac_outputs.mat');
if ~exist(cachePath, 'file')
    error('run_beamforming_sweep:cacheMissing', ...
        ['LAC output cache not found at:\n  %s\nPlease run run_paper_figures.m first.\n' ...
         '未找到 LAC 输出缓存，请先运行 run_paper_figures.m。'], cachePath);
end
fprintf('Loading cached LAC outputs / 加载 LAC 输出缓存:\n  %s\n', cachePath);
cache    = load(cachePath);
srir_orig = cache.srir_orig;
srir_lac  = cache.srir_new;
srir_tar  = cache.srir_tar;
fs        = cache.fs;

% --- 2) build decoder once (paper script reuses what's in base) -----------
build_ambisonic_decoder(repoRoot);

% Push LAC inputs + silence flags to base (kept for every sweep iteration).
assignin('base', 'idx_LS_src_orig',      cache.idx_LS_src_orig);
assignin('base', 'idx_LS_rec_orig',      cache.idx_LS_rec_orig);
assignin('base', 'idx_LS_src_tar',       cache.idx_LS_src_tar);
assignin('base', 'idx_LS_rec_tar',       cache.idx_LS_rec_tar);
assignin('base', 'flag_use_beamforming', true);
assignin('base', 'save_figs',            false);
assignin('base', 'flag_render_binaural', false);
assignin('base', 'repoRoot',             repoRoot);

% Hide figures the paper script creates internally during each iteration.
prevVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
restoreVis = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));   %#ok<NASGU>

% --- 3) LAC baseline losses (one-off) -------------------------------------
fprintf('\nComputing LAC baseline losses / 计算 LAC 基线损失 ...\n');
[L_lac, P_lac] = compute_losses(srir_lac, srir_tar, fs);
fprintf('  LAC: L_time=%.4f, L_spec=%.3f dB, L_pwd=%.4f\n', ...
    L_lac.time, L_lac.spec, L_lac.pwd);

% --- 4) Coarse grid sweep --------------------------------------------------
% Coarse grid extended toward small beta because the first sweep showed
% the optimum at the lower beta edge (best was beta=0.5).
% 把 grid 往小 beta 方向延伸（首轮 best 在 beta=0.5 边界）。
betas  = [0.2  0.3  0.4  0.5  0.7];
sigmas = [0.8  1.2  1.8  2.5  3.5];
nB = numel(betas); nS = numel(sigmas);

fprintf('\nSweeping over %d x %d grid (%d evaluations) ...\n', nB, nS, nB*nS);
srir_grid    = cell(nB, nS);
L_time_grid  = nan(nB, nS);
L_spec_grid  = nan(nB, nS);
L_pwd_grid   = nan(nB, nS);
L_comp_grid  = nan(nB, nS);     % composite (normalised against LAC)
P_grid       = cell(nB, nS);
gd_ref       = [];

t0 = tic;
for ib = 1:nB
    for is = 1:nS
        fprintf('  [%2d/%2d]  beta=%.2f  sigma=%.2f ... ', ...
            (ib-1)*nS + is, nB*nS, betas(ib), sigmas(is));
        tA = tic;
        srir_beam = run_lac_with_beam(betas(ib), sigmas(is));
        [Lk, Pk, gd_ref] = compute_losses(srir_beam, srir_tar, fs);
        srir_grid{ib, is}   = srir_beam;
        P_grid  {ib, is}    = Pk;
        L_time_grid(ib, is) = Lk.time;
        L_spec_grid(ib, is) = Lk.spec;
        L_pwd_grid (ib, is) = Lk.pwd;
        L_comp_grid(ib, is) = mean([Lk.time / L_lac.time, ...
                                    Lk.spec / L_lac.spec, ...
                                    Lk.pwd  / L_lac.pwd]);
        fprintf('comp=%.3f  (%.1fs)\n', L_comp_grid(ib, is), toc(tA));
    end
end
fprintf('Grid sweep done in %.1f s.\n', toc(t0));

% --- 5) Pick best on grid; decide whether refinement is needed ------------
[L_best_grid, idx_best] = min(L_comp_grid(:));
[ib_best, is_best]      = ind2sub([nB nS], idx_best);
beta_best  = betas(ib_best);
sigma_best = sigmas(is_best);
fprintf(['\nGrid best:  beta=%.2f, sigma=%.2f\n' ...
         '  composite loss = %.3f (LAC normalised; <1.0 beats LAC)\n'], ...
         beta_best, sigma_best, L_best_grid);

% Always run fminsearch refinement starting from the grid best, even if
% the grid result already beats LAC. This catches optima that fall between
% grid nodes or near the edge of the swept range.
% 始终从 grid best 出发跑一次 fminsearch 精化 —— 不管 grid best 多好，
% 因为真实最优可能落在 grid 节点之间或边界附近。
do_refine = true;
fprintf('Always running fminsearch refinement from grid best ...\n');

% --- 6) Optional fminsearch refinement ------------------------------------
% Variable transform to keep beta in (0,1) and sigma > 0 without bounds:
%   beta  = 1/(1+exp(-x1))     (sigmoid)
%   sigma = exp(x2)            (positive)
% fminsearch then runs in unbounded R^2.
% 变量变换：beta 通过 sigmoid 限制在 (0,1)，sigma 通过 exp 保正；fminsearch
% 在无界 R^2 上跑，无需 toolbox。
if do_refine
    x0   = [log(beta_best / (1 - beta_best)), log(sigma_best)];
    fobj = @(x) refine_objective(x, srir_tar, fs, L_lac);
    opts = optimset('Display','iter', 'MaxFunEvals', 25, 'MaxIter', 25, ...
                    'TolX', 1e-3, 'TolFun', 1e-3);
    [x_opt, L_final] = fminsearch(fobj, x0, opts);
    beta_final  = 1 / (1 + exp(-x_opt(1)));
    sigma_final = exp(x_opt(2));
    fprintf('fminsearch final: beta=%.3f, sigma=%.3f, comp=%.3f\n', ...
        beta_final, sigma_final, L_final);
else
    beta_final  = beta_best;
    sigma_final = sigma_best;
    L_final     = L_best_grid;
end

% --- 7) Compute the final SRIR with the chosen (beta, sigma) --------------
fprintf('\nComputing final SRIR at chosen (beta=%.3f, sigma=%.3f) ...\n', ...
    beta_final, sigma_final);
srir_final = run_lac_with_beam(beta_final, sigma_final);
[L_final_pieces, P_final] = compute_losses(srir_final, srir_tar, fs);

improvement = struct( ...
    'time', 100*(L_lac.time - L_final_pieces.time)/L_lac.time, ...
    'spec', 100*(L_lac.spec - L_final_pieces.spec)/L_lac.spec, ...
    'pwd',  100*(L_lac.pwd  - L_final_pieces.pwd )/L_lac.pwd);

% --- 8) Plot 5x5 grid of spherical power maps -----------------------------
fprintf('\nPlotting sweep grid / 绘制 sweep 大宫格 ...\n');
figDir = fullfile(repoRoot, 'figures');
if ~isfolder(figDir); mkdir(figDir); end

figh1 = figure('Color','w','Position',[40 40 1700 1200],'Visible','on');
for ib = 1:nB
    for is = 1:nS
        subplot(nB, nS, (ib-1)*nS + is);
        heatmap_plot(rad2deg(gd_ref(:,1)), rad2deg(gd_ref(:,2)), P_grid{ib,is});
        colormap(flipud(bone)); clim([0 3]);
        set(gca, 'FontSize', 8); xlabel(''); ylabel('');
        winnerMark = '';
        if ib == ib_best && is == is_best; winnerMark = '  *'; end
        title(sprintf('\\beta=%.2g, \\sigma=%.2g (L=%.2f)%s', ...
            betas(ib), sigmas(is), L_comp_grid(ib,is), winnerMark), 'FontSize', 9);
    end
end
sgtitle(sprintf(['Beamforming sweep over (\\beta, \\sigma).  ' ...
    'Composite loss vs LAC: <1.0 better, >1.0 worse.   ' ...
    'Grid best * = (\\beta=%.2g, \\sigma=%.2g, L=%.2f)'], ...
    beta_best, sigma_best, L_best_grid), 'FontWeight', 'bold');
gridPdf = fullfile(figDir, 'beamforming_sweep_grid.pdf');
gridPng = fullfile(figDir, 'beamforming_sweep_grid.png');
exportgraphics(figh1, gridPdf);
exportgraphics(figh1, gridPng, 'Resolution', 180);

% --- 9) Final 4-panel comparison ------------------------------------------
figh2 = figure('Color','w','Position',[60 60 1900 380],'Visible','on');

subplot(1,4,1);
P = compute_pwd(srir_orig, fs);
heatmap_plot(rad2deg(gd_ref(:,1)), rad2deg(gd_ref(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca,'FontSize',11);
title('Original SRIR');

subplot(1,4,2);
heatmap_plot(rad2deg(gd_ref(:,1)), rad2deg(gd_ref(:,2)), P_lac);
colormap(flipud(bone)); clim([0 3]); set(gca,'FontSize',11);
title(sprintf('LAC extrapolation  (L_{time}=%.3f, L_{spec}=%.2fdB, L_{pwd}=%.3f)', ...
    L_lac.time, L_lac.spec, L_lac.pwd));

subplot(1,4,3);
heatmap_plot(rad2deg(gd_ref(:,1)), rad2deg(gd_ref(:,2)), P_final);
colormap(flipud(bone)); clim([0 3]); set(gca,'FontSize',11);
title(sprintf('Beamforming-refined  \\beta=%.2f, \\sigma=%.2f  (comp=%.3f)', ...
    beta_final, sigma_final, L_final));

subplot(1,4,4);
P = compute_pwd(srir_tar, fs);
heatmap_plot(rad2deg(gd_ref(:,1)), rad2deg(gd_ref(:,2)), P);
colormap(flipud(bone)); clim([0 3]); set(gca,'FontSize',11);
title('Target (ground truth)');

cb = colorbar('Position', [0.945 0.18 0.010 0.7]);
ylabel(cb, 'Normalised power (dB)');

finalPdf = fullfile(figDir, 'beamforming_final_4panel.pdf');
finalPng = fullfile(figDir, 'beamforming_final_4panel.png');
exportgraphics(figh2, finalPdf);
exportgraphics(figh2, finalPng, 'Resolution', 200);

% --- 10) CSV with all losses ----------------------------------------------
csvPath = fullfile(figDir, 'beamforming_sweep_losses.csv');
fid = fopen(csvPath, 'w');
fprintf(fid, 'beta,sigma,L_time,L_spec_dB,L_pwd,L_composite_vs_LAC\n');
for ib = 1:nB
    for is = 1:nS
        fprintf(fid, '%.3f,%.3f,%.6f,%.6f,%.6f,%.6f\n', ...
            betas(ib), sigmas(is), ...
            L_time_grid(ib,is), L_spec_grid(ib,is), ...
            L_pwd_grid(ib,is),  L_comp_grid(ib,is));
    end
end
fprintf(fid, 'FINAL,FINAL,%.6f,%.6f,%.6f,%.6f\n', ...
    L_final_pieces.time, L_final_pieces.spec, L_final_pieces.pwd, L_final);
fclose(fid);

% --- 11) Final report ------------------------------------------------------
fprintf('\n================ BEAMFORMING SWEEP SUMMARY ================\n');
fprintf('LAC baseline:  L_time=%.4f   L_spec=%.3f dB   L_pwd=%.4f\n', ...
    L_lac.time, L_lac.spec, L_lac.pwd);
fprintf('FINAL (beta=%.3f, sigma=%.3f):\n', beta_final, sigma_final);
fprintf('               L_time=%.4f   L_spec=%.3f dB   L_pwd=%.4f\n', ...
    L_final_pieces.time, L_final_pieces.spec, L_final_pieces.pwd);
fprintf('Improvement over LAC:\n');
fprintf('               time : %+.1f%%   (positive = better)\n', improvement.time);
fprintf('               spec : %+.1f%%\n',                       improvement.spec);
fprintf('               pwd  : %+.1f%%\n',                       improvement.pwd);
fprintf('Composite loss vs LAC: %.3f (1.0 = same; <1.0 = better)\n', L_final);
fprintf('\nFiles written:\n');
fprintf('  %s\n', gridPdf);
fprintf('  %s\n', finalPdf);
fprintf('  %s\n', csvPath);
fprintf('===========================================================\n');


% =========================================================================
% local helpers / 局部函数
% =========================================================================

function v = refine_objective(x, srir_tar, fs, L_lac)
% fminsearch objective: x = [logit(beta), log(sigma)]; returns composite
% loss normalised to LAC baseline (lower = better).
% fminsearch 目标函数：x = [logit(beta), log(sigma)]，返回归一化到 LAC
% 基线的复合损失（越小越好）。
b  = 1 / (1 + exp(-x(1)));
s  = exp(x(2));
sr = run_lac_with_beam(b, s);
Lh = compute_losses(sr, srir_tar, fs);
v  = mean([Lh.time / L_lac.time, Lh.spec / L_lac.spec, Lh.pwd / L_lac.pwd]);
fprintf('    refine eval  beta=%.3f  sigma=%.3f  comp=%.3f\n', b, s, v);
end


function srir_new = run_lac_with_beam(beta_val, sigma_val)
% Run the paper LAC script once with the given (beta, sigma), silently,
% and return the resulting srir_new from the base workspace.
% 用给定 (beta, sigma) 静默跑一次论文 LAC 脚本，返回 srir_new。
assignin('base', 'beam_beta',  beta_val);
assignin('base', 'beam_sigma', sigma_val);
evalin('base', 'srir_extrap_imgSrc_subspace_LAC');
close all hidden;
srir_new = evalin('base', 'srir_new');
end


function [L, P_pwd, grid_dirs] = compute_losses(srir, srir_tar, fs)
% Three-component loss against the target SRIR. Returns the PWD power map
% as a by-product so the caller does not need to recompute for plotting.
% 计算三项损失（vs target），顺便返回 PWD 功率图供绘图复用。
[P_pwd, grid_dirs] = compute_pwd(srir,     fs);
[P_tar, ~        ] = compute_pwd(srir_tar, fs);

% Time-domain NRMSE on W channel (channel 1 in ACN order).
w  = srir(:,1);  wt = srir_tar(:,1);
n  = min(numel(w), numel(wt));
L.time = norm(w(1:n) - wt(1:n)) / max(norm(wt(1:n)), eps);

% 1/3-oct W-channel log-spectrum mean abs deviation.
L.spec = third_oct_spec_diff_dB(w(1:n), wt(1:n), fs);

% Spherical PWD L1 distance after L1-normalisation (probability-like).
pb = P_pwd / max(sum(P_pwd), eps);
pt = P_tar / max(sum(P_tar), eps);
L.pwd = sum(abs(pb - pt));
end


function [P_pwd, grid_dirs] = compute_pwd(srir, fs)
% PWD power map identical in spirit to the paper's get_pwd (Fig. 6).
% 与论文 get_pwd 等价的 PWD 功率图。
degRes  = 2;
order   = sqrt(size(srir, 2)) - 1;
nSrc    = 7;
nSmp    = round(fs * 0.3);
fHp     = 3000;
kappa   = 40;
grid_dirs = grid2dirs(degRes, degRes, 0, 0);
P_src     = diag(ones(nSmp, 1));
[~, filtHi, ~] = ambisonic_crossover(fHp, fs);
Y_src     = filter(filtHi, 1, srir(1:nSmp, :));
stVec     = Y_src';
sphCOV    = stVec * P_src * stVec' + eye((order + 1)^2) / (4 * pi);
[P_pwd, ~, ~] = sphPWDmap_compat(sphCOV, grid_dirs, nSrc, kappa);
end


function L = third_oct_spec_diff_dB(x, y, fs)
% Average absolute log-spectrum difference (dB) integrated over 1/3-oct
% bands between 100 Hz and ~16 kHz.
% 100 Hz–16 kHz 区间内，按 1/3 oct 频带积分的对数频谱差异（dB）。
N    = max(numel(x), numel(y));
Nfft = 2^nextpow2(N);
X    = abs(fft(x, Nfft));   X = X(1:Nfft/2);
Y    = abs(fft(y, Nfft));   Y = Y(1:Nfft/2);
freqs = (0:Nfft/2-1).' * (fs / Nfft);

% 1/3 oct centre frequencies
fc = 100 * 2.^((0:1/3:7)');
fc = fc(fc < fs/2 * 0.95);

L_acc = 0; nB = 0;
for i = 1:numel(fc)
    fl = fc(i) / 2^(1/6);
    fu = fc(i) * 2^(1/6);
    sel = freqs >= fl & freqs <= fu;
    if any(sel)
        ex = mean(X(sel).^2);
        ey = mean(Y(sel).^2);
        L_acc = L_acc + abs(10*log10(max(ex,eps)) - 10*log10(max(ey,eps)));
        nB    = nB + 1;
    end
end
L = L_acc / max(nB, 1);
end
