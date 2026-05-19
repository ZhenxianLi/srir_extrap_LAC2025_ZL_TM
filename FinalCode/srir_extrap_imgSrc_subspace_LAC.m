%% SRIR_EXTRAP_IMGSRC_SUBSPACE_LAC
%  Spatial room impulse response (SRIR) extrapolation using image-source
%  modelling and direct/residual subspace decomposition.
%  基于镜像源法 + 直达/残差子空间分解的空间房间脉冲响应外推。
%
%  This script is the one referenced in the LAC-25 paper. It is intended to
%  be invoked from `run_paper_figures.m` (which sets up the MATLAB path,
%  downloads the dataset, loads the Ambisonic decoder, and chdir-s into
%  the `figures/` folder before running this).
%  本脚本由 run_paper_figures.m 调用：后者负责加 path、下载数据、加载
%  Ambisonic 解码器，并把工作目录切到 figures/。
%
%  Paper / 论文:
%    Zhenxian Li and Thomas McKenzie,
%    "Spatial Room Impulse Response Extrapolation using the Image Source Method,"
%    Proc. 19th Linux Audio Conference (LAC-25), 2025.
%
%  Author / 作者: Zhenxian Li
% =========================================================================

% NOTE / 注意:
%   We deliberately do NOT call `clear all` / `clc` / `close all` here —
%   the launcher already set up the workspace and prepared the decoder.
%   This script can also be reused by other launchers (e.g. beamforming).
%   故意不在此处 clear/clc/close —— 上层启动脚本已处理。

%% Settings / 参数
% Honour any value pre-set by an outer launcher (run_beamforming.m sets
% save_figs=false to avoid overwriting paper PDFs during its second pass).
% 尊重外层 launcher 预设的值（run_beamforming.m 会预先把 save_figs 设为
% false，避免第二次跑时覆盖论文 PDF）。
if ~exist('save_figs', 'var'); save_figs = true; end
dim_room               = [7.87 5.75 2.91];   % from the ArXiv paper / 自论文
ord_imgsrc             = 2;                  % image-source order / 镜像源阶数

% Original (measured) and target (extrapolated) source/receiver indices.
% 原始（实测）和目标（外推）声源/接收点编号。论文示例：S2-R4 -> S3-R2
idx_LS_src_orig = 2;
idx_LS_rec_orig = 4;
idx_LS_src_tar  = 3;
idx_LS_rec_tar  = 2;

% Arrival-extraction windows (samples / ms). / 到达分量提取窗口长度
delay_arrival_initial_samp     = 10;
arrival_dur_ms                 = 1;
delay_arrival_initial_samp_DS  = 50;
arrival_dur_ms_DS              = 4;

% Behaviour flags / 行为开关
% Same "honour external value" pattern — run_beamforming.m sets
% flag_render_binaural=false so it does not overwrite the LAC audio with
% the beamforming-refined audio.
% 同样的"尊重外部值"模式：run_beamforming.m 设 flag_render_binaural=false
% 以免用波束细化版的音频覆盖 LAC 版的音频。
if ~exist('flag_subspace_decomp', 'var'); flag_subspace_decomp = true; end
if ~exist('flag_render_binaural', 'var'); flag_render_binaural = true; end   % false 可跳过音频导出

% Beamforming-refined extension (set by run_beamforming.m before re-invoking
% this script; default is the plain LAC paper method).
% beamforming-refined 扩展：默认 false 走论文 LAC；run_beamforming.m 会在
% base workspace 里把 flag 置为 true 再次调用本脚本，得到 beamforming 版本。
if ~exist('flag_use_beamforming', 'var'); flag_use_beamforming = false; end
if ~exist('beam_beta',  'var');           beam_beta  = 0.8;             end
if ~exist('beam_sigma', 'var');           beam_sigma = 1.2;             end

%% Locate dataset / 定位数据集
% repoRoot is set by run_paper_figures.m via setupPaths; fall back if missing
if ~exist('repoRoot', 'var') || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
irPath1 = fullfile(repoRoot, '6dof_SRIRs_eigenmike_SH');
irName1 = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';
sofaFull = fullfile(irPath1, irName1);
assert(exist(sofaFull, 'file') == 2, ...
    'SOFA file missing: %s — run downloadSrirDataset first.\n找不到 SOFA 数据，请先运行 downloadSrirDataset。', sofaFull);

sofa1 = SOFAload(sofaFull);
fs    = sofa1.Data.SamplingRate;

%% Get original and target SRIRs / 取出原始与目标 SRIR
idx_LS_srcrec_orig = idx_LS_rec_orig*3 - 3 + idx_LS_src_orig;  % flat index
idx_LS_srcrec_tar  = idx_LS_rec_tar *3 - 3 + idx_LS_src_tar;

srir_orig = squeeze(sofa1.Data.IR(idx_LS_srcrec_orig, :, :))';
srir_tar  = squeeze(sofa1.Data.IR(idx_LS_srcrec_tar,  :, :))';

pos_src_orig  = sofa1.SourcePosition(idx_LS_srcrec_orig, :);
pos_list_orig = sofa1.ListenerPosition(idx_LS_srcrec_orig, :);
pos_src_tar   = sofa1.SourcePosition(idx_LS_srcrec_tar,  :);
pos_list_tar  = sofa1.ListenerPosition(idx_LS_srcrec_tar,  :);

% Plot geometry / 绘制几何
% Tricky: SOFAplotGeometry's switch statement *does* handle SingleRoomDRIR
% (line ~107 in SOFAtoolbox 2.4), but it crashes inside because the new
% SOFAtoolbox no longer ships the SingleRoomDRIR.csv convention file
% (the convention was renamed to SingleRoomSRIR upstream). The plotting
% branch for SingleRoomSRIR is identical, and our SingleRoomDRIR fields
% are a subset of SingleRoomSRIR, so swap the convention name on a shim
% copy and re-try. Final fallback is a manual scatter + view-vector plot.
% 玄机：SOFAplotGeometry 的 switch 里其实包含 SingleRoomDRIR 分支（约 107 行），
% 但新版 SOFAtoolbox 2.x 不再附带 SingleRoomDRIR.csv，所以内部 SOFAexpand
% 取 convention 时会崩。该 convention 在上游被改名为 SingleRoomSRIR；其
% 绘图逻辑分支完全一样，字段也兼容。所以这里在副本里把 convention 改名后
% 再调一次，仍失败则退回到手动散点 + 朝向箭头。
plotted = false;
try
    SOFAplotGeometry(sofa1);
    plotted = true;
catch
    sofaShim = sofa1;
    if strcmp(sofaShim.GLOBAL_SOFAConventions, 'SingleRoomDRIR')
        sofaShim.GLOBAL_SOFAConventions = 'SingleRoomSRIR';
        sofaShim.GLOBAL_SOFAConventionsVersion = '1.0';
        try
            SOFAplotGeometry(sofaShim);
            plotted = true;
        catch
            % fall through to manual fallback below / 落到下面手动 fallback
        end
    end
end
if ~plotted
    warning('Geometry:fallback', ...
        ['SOFAplotGeometry could not render; drawing scatter + view arrows.\n' ...
         'SOFAplotGeometry 无法绘制，已退回到手动散点 + 朝向箭头。']);
    figure; hold on; grid on; axis equal;
    plot3(sofa1.SourcePosition(:,1),   sofa1.SourcePosition(:,2),   sofa1.SourcePosition(:,3),   'b.', 'MarkerSize', 10);
    plot3(sofa1.ListenerPosition(:,1), sofa1.ListenerPosition(:,2), sofa1.ListenerPosition(:,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    if isfield(sofa1, 'ListenerView')
        LP = sofa1.ListenerPosition;  LV = sofa1.ListenerView;
        if size(LV,1) == 1; LV = repmat(LV, size(LP,1), 1); end
        quiver3(LP(:,1), LP(:,2), LP(:,3), LV(:,1), LV(:,2), LV(:,3), 0.3, 'Color', 'r');
    end
    if isfield(sofa1, 'SourceView')
        SP = sofa1.SourcePosition;  SV = sofa1.SourceView;
        if size(SV,1) == 1; SV = repmat(SV, size(SP,1), 1); end
        quiver3(SP(:,1), SP(:,2), SP(:,3), SV(:,1), SV(:,2), SV(:,3), 0.3, 'Color', 'b');
    end
    xlabel('x / m'); ylabel('y / m'); zlabel('z / m');
    legend({'SourcePosition', 'ListenerPosition', 'ListenerView', 'SourceView'}, 'Location', 'best');
    view(2);
end
line([pos_src_orig(1) pos_list_orig(1)], [pos_src_orig(2) pos_list_orig(2)], 'Color','blue');
line([pos_src_tar(1)  pos_list_tar(1) ], [pos_src_tar(2)  pos_list_tar(2) ], 'Color','red');
plot_time_vector = 0:1/fs:(size(srir_orig,1)/fs - 1/fs);

%% Intensity-vector DoA of measured direct sound / 实测直达声 DoA
[delay_DS_orig_meas, ~, ~, ~] = findDirectSound(srir_orig(:,1));
ds_start_samp = delay_DS_orig_meas - 0.001*fs;
ds_end_samp   = ds_start_samp + 0.004*fs;
dir_DS_meas_orig = doa_iv(srir_orig(ds_start_samp:ds_end_samp, :));

%% Subspace decomposition / 子空间分解
if flag_subspace_decomp
    blockLenSmp = 32;  hopSizeSmp = blockLenSmp/8;  kappa = 2.5;
    numBlocksGsvSumAvg       = 32;
    residualEstimateLengthMs = 20;
    decompositionTimeLimitMs = 100;
    numBlocksSmoothThresh    = 1;
    [srir_direct, srir_resid, ~, ~, ~, ~, ~] = ...
        srirSubspaceDecomp(srir_orig, fs, blockLenSmp, hopSizeSmp, kappa, ...
            numBlocksGsvSumAvg, residualEstimateLengthMs, ...
            decompositionTimeLimitMs, numBlocksSmoothThresh);
end

%% Image-source detection (original) / 镜像源计算（原始位置）
[dist_imgsrc_orig, delay_imgsrc_orig, dir_imgsrc_orig, ...
 dist_DS_imgsrc_orig, delay_DS_imgsrc_orig, dir_DS_imgsrc_orig] = ...
    imgsrc_calculate(pos_list_orig, pos_src_orig, ord_imgsrc, dim_room, fs);
view([0 90]);
xlim([-15 15 + dim_room(1)]);  ylim([-12 12 + dim_room(2)]);
set(gca, 'FontSize', 12);
pbaspect([15 + dim_room(1) + 15, 12 + dim_room(2) + 12, 1]);
box on;
if save_figs; exportgraphics(gcf, 'srir_extrap_imgsrc_orig.pdf'); end

%% Image-source detection (target) / 镜像源计算（目标位置）
[dist_imgsrc_tar, delay_imgsrc_tar, dir_imgsrc_tar, ...
 dist_DS_imgsrc_tar, delay_DS_imgsrc_tar, dir_DS_imgsrc_tar] = ...
    imgsrc_calculate(pos_list_tar, pos_src_tar, ord_imgsrc, dim_room, fs);
view([0 90]);
xlim([-15 15 + dim_room(1)]);  ylim([-12 12 + dim_room(2)]);
set(gca, 'FontSize', 12);
pbaspect([15 + dim_room(1) + 15, 12 + dim_room(2) + 12, 1]);
box on;
if save_figs; exportgraphics(gcf, 'srir_extrap_imgsrc_tar.pdf'); end

% Sort arrivals by mean distance / 按平均距离排序，最重要的先处理
[~, idx] = sort(mean([dist_imgsrc_orig; dist_imgsrc_tar]));
if ord_imgsrc > 1  % drop the duplicate direct-sound path at order >=2
    idx = idx(2:end);
end
dist_imgsrc_orig  = dist_imgsrc_orig(idx);
delay_imgsrc_orig = delay_imgsrc_orig(idx);
dir_imgsrc_orig   = dir_imgsrc_orig(:, idx);
dist_imgsrc_tar   = dist_imgsrc_tar(idx);
delay_imgsrc_tar  = delay_imgsrc_tar(idx);
dir_imgsrc_tar    = dir_imgsrc_tar(:, idx);

%% Align measured DS direction with image-source DS direction
%  对齐实测直达声方向与镜像源直达声方向
dir_diff_DS = rad2deg(angdiff(deg2rad(dir_DS_meas_orig), deg2rad(dir_DS_imgsrc_orig)));
dir_diff_threshold = 30;  % deg
if abs(dir_diff_DS(1,1)) > dir_diff_threshold
    dir_diff_DS_rnd = dir_diff_threshold * round(dir_diff_DS / dir_diff_threshold);
    fprintf(['Direction mismatch %g° -> counter-rotating %g°\n' ...
             '方向差 %g°，反向旋转 %g° 校正\n'], ...
            dir_diff_DS(1,1), dir_diff_DS_rnd(1,1), dir_diff_DS(1,1), dir_diff_DS_rnd(1,1));
    dir_DS_imgsrc_orig = dir_DS_imgsrc_orig - dir_diff_DS_rnd;
    dir_DS_imgsrc_tar  = dir_DS_imgsrc_tar  - dir_diff_DS_rnd;
    dir_imgsrc_orig    = dir_imgsrc_orig    - dir_diff_DS_rnd;
    dir_imgsrc_tar     = dir_imgsrc_tar     - dir_diff_DS_rnd;
end

%% Align measured DS time with image-source DS time / 对齐到达时间
delay_diff_DS = delay_DS_orig_meas - delay_DS_imgsrc_orig;
delay_diff_threshold_samp = 0.1/1000 * fs;
if abs(delay_diff_DS) > delay_diff_threshold_samp
    fprintf(['Time mismatch %.3f ms -> counter-delaying %d samples\n' ...
             '到达时间差 %.3f ms，补偿 %d 个采样\n'], ...
            delay_diff_DS/fs*1000, round(delay_diff_DS), ...
            delay_diff_DS/fs*1000, round(delay_diff_DS));
    delay_DS_imgsrc_orig = delay_DS_imgsrc_orig + round(delay_diff_DS);
    delay_imgsrc_orig    = delay_imgsrc_orig    + round(delay_diff_DS);
    delay_DS_imgsrc_tar  = delay_DS_imgsrc_tar  + round(delay_diff_DS);
    delay_imgsrc_tar     = delay_imgsrc_tar     + round(delay_diff_DS);
end

%% Transform parameters (orig -> tar) / 计算原始到目标的变换参数
dir_diff_DS_transform = rad2deg(angdiff(deg2rad(dir_DS_imgsrc_orig), deg2rad(dir_DS_imgsrc_tar)));
dir_diff_transform    = rad2deg(angdiff(deg2rad(dir_imgsrc_orig),    deg2rad(dir_imgsrc_tar)));
delay_diff_DS_transform = round(delay_DS_imgsrc_tar - delay_DS_imgsrc_orig);
delay_diff_transform    = round(delay_imgsrc_tar    - delay_imgsrc_orig);
dist_diff_DS_transform  = dist_DS_imgsrc_orig ./ dist_DS_imgsrc_tar;
dist_diff_transform     = dist_imgsrc_orig    ./ dist_imgsrc_tar;

plotcutoff = round(fs/20);

%% Extrapolation: window-out, rotate, delay, gain, re-insert
%  外推主循环：截窗 → 旋转 → 延时 → 增益 → 重新插入
gainpowerDS = 1;  gainpower = 1;
if flag_subspace_decomp
    srir_arrival_removed = srir_direct;
    srir_unAlt           = srir_direct;
else
    srir_arrival_removed = srir_orig;
    srir_unAlt           = srir_orig;
end

% Early reflections / 早期反射
for i = 1:length(delay_imgsrc_orig)
    [srir_arrival, srir_arrival_removed] = srir_arrival_remove( ...
        srir_arrival_removed, srir_unAlt, fs, ...
        delay_imgsrc_orig(i) - delay_arrival_initial_samp, arrival_dur_ms);

    dir_arr_orig = dir_imgsrc_orig(:, i);
    dir_arr_tar  = dir_imgsrc_tar(:,  i);

    % Beamforming-refined: filter the windowed arrival in the direction of
    % the image-source so that energy from concurrent reflections inside
    % the same time window is suppressed before rotation+re-insertion.
    % 波束细化：用镜像源方向对该窗做方向选择性过滤，先剔除窗内其它方向
    % 反射的干扰，再做旋转和重新插入。
    if flag_use_beamforming
        srir_arrival = beam_filter(srir_arrival, ...
            dir_arr_orig(1), dir_arr_orig(2), beam_beta, beam_sigma);
    end

    dir_rot = rad2deg(angdiff(deg2rad(dir_arr_orig), deg2rad(dir_arr_tar)));
    srir_arrival_rot = rotateHOA_N3D(srir_arrival, dir_rot(1), dir_rot(2), 0);

    s = delay_imgsrc_tar(i) - delay_arrival_initial_samp + 1;
    e = s + size(srir_arrival, 1) - 1;
    srir_arrival_removed(s:e, :) = srir_arrival_removed(s:e, :) + ...
        srir_arrival_rot .* dist_diff_transform(i) .^ gainpower;
end

% Direct sound / 直达声
[srir_arrival, srir_arrival_removed] = srir_arrival_remove( ...
    srir_arrival_removed, srir_unAlt, fs, ...
    delay_DS_imgsrc_orig - delay_arrival_initial_samp_DS, arrival_dur_ms_DS);

% Beamforming-refined: same direction-selective filter on the DS window.
% 波束细化：对直达声窗也做一次方向过滤。
if flag_use_beamforming
    srir_arrival = beam_filter(srir_arrival, ...
        dir_DS_imgsrc_orig(1), dir_DS_imgsrc_orig(2), beam_beta, beam_sigma);
end

dir_rot = rad2deg(angdiff(deg2rad(dir_DS_imgsrc_orig), deg2rad(dir_DS_imgsrc_tar)));
srir_arrival_rot = rotateHOA_N3D(srir_arrival, dir_rot(1), dir_rot(2), 0);

srir_arrival_removed(1:delay_DS_imgsrc_tar - delay_arrival_initial_samp_DS + 1, :) = 0;
srir_resid          (1:delay_DS_imgsrc_tar - delay_arrival_initial_samp_DS + 1, :) = 0;

s = delay_DS_imgsrc_tar - delay_arrival_initial_samp_DS + 1;
e = s + size(srir_arrival, 1) - 1;
srir_arrival_removed(s:e, :) = srir_arrival_removed(s:e, :) + ...
    srir_arrival_rot .* dist_diff_DS_transform .^ gainpowerDS;

% Combine with diffuse tail / 与扩散残差合并
if flag_subspace_decomp
    srir_new = srir_arrival_removed + srir_resid;
else
    srir_new = srir_arrival_removed;
end

%% ===================== EVALUATION / 评价 =================================

% ---- Time-domain comparison (Fig. 3) / 时域对比 ---------------------------
figure;
subplot(3,1,1);
plot(plot_time_vector(1:plotcutoff), db(abs(srir_orig(1:plotcutoff,1))));
title(sprintf('Original. RMS = %.4g dB', db(rms(srir_orig(:,1)))));
ylim([-40 0]); ylabel('Magnitude (dB)'); xlabel('Time (s)');

subplot(3,1,2);
plot(plot_time_vector(1:plotcutoff), db(abs(srir_new(1:plotcutoff,1))));
title(sprintf('Extrapolated. RMS = %.4g dB', db(rms(srir_new(:,1)))));
ylim([-40 0]); ylabel('Magnitude (dB)'); xlabel('Time (s)');

subplot(3,1,3);
plot(plot_time_vector(1:plotcutoff), db(abs(srir_tar(1:plotcutoff,1))));
title(sprintf('Target. RMS = %.4g dB', db(rms(srir_tar(:,1)))));
ylim([-40 0]); ylabel('Magnitude (dB)'); xlabel('Time (s)');

% ---- Comparison overlays (Fig. 3 panels in paper) -----------------------
figure;
plot(plot_time_vector(1:plotcutoff), db(abs(srir_orig(1:plotcutoff,1)))); hold on;
plot(plot_time_vector(1:plotcutoff), db(abs(srir_new (1:plotcutoff,1))));
legend({'Original','Extrapolated'}); ylim([-30 0]);
ylabel('Magnitude (dB)'); xlabel('Time (s)');
pbaspect([3 1 1]); grid on; box on;
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_comparison_orig.pdf'); end

figure;
plot(plot_time_vector(1:plotcutoff), db(abs(srir_tar(1:plotcutoff,1)))); hold on;
plot(plot_time_vector(1:plotcutoff), db(abs(srir_new(1:plotcutoff,1))));
legend({'Target','Extrapolated'}); ylim([-30 0]);
ylabel('Magnitude (dB)'); xlabel('Time (s)');
pbaspect([3 1 1]); grid on; box on;
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_comparison_extrap.pdf'); end

%% Binaural render / 双耳渲染
% SH_ambisonic_binaural_decoder is expected to live in the base workspace
% (loaded by run_paper_figures.m via load_ambisonic_configuration).
% SH_ambisonic_binaural_decoder 由 run_paper_figures.m 加载到 base 工作区。
if ~exist('SH_ambisonic_binaural_decoder', 'var')
    SH_ambisonic_binaural_decoder = evalin('base', 'SH_ambisonic_binaural_decoder');
end

brir_orig = binSound(srir_orig, SH_ambisonic_binaural_decoder);
brir_new  = binSound(srir_new,  SH_ambisonic_binaural_decoder);
brir_tar  = binSound(srir_tar,  SH_ambisonic_binaural_decoder);

%% Predicted Binaural Colouration (PBC-2) / 预测双耳音色 (PBC-2)
% Needs `mckenzie2025` — the PBC-2 model used in the LAC-25 paper. As of
% AMT 1.6.0 it is NOT yet merged upstream (only mckenzie2021/2022 are).
% To reproduce this block, obtain mckenzie2025.m from Thomas McKenzie and
% drop it into external/amtoolbox/<...>/models/. Otherwise this block is
% skipped with a warning — everything else (figures + audio) still runs.
% 需要 `mckenzie2025`（即论文里 PBC-2 模型）。截至 AMT 1.6.0 该模型尚未
% 合入主线（只有 mckenzie2021/2022）。要复现此段，请向 Thomas McKenzie
% 索取 mckenzie2025.m 并放入 external/amtoolbox/<...>/models/。否则
% 此段会被跳过并给出警告，不影响其它图与音频生成。
if exist('mckenzie2025', 'file') == 2
    brir_o = reshape(brir_orig, [], 1, 2);
    brir_n = reshape(brir_new,  [], 1, 2);
    brir_t = reshape(brir_tar,  [], 1, 2);
    settings.smGL1 = 0;
    pbc2_orig = mckenzie2025(brir_o, brir_t, settings);
    pbc2_new  = mckenzie2025(brir_n, brir_t, settings);
    fprintf('PBC (orig -> tar) = %.3g\n',  pbc2_orig);
    fprintf('PBC (new  -> tar) = %.3g\n',  pbc2_new);
else
    warning('mckenzie2025:notFound', ...
        ['mckenzie2025 not on MATLAB path — PBC-2 colouration block skipped.\n' ...
         'Note: as of AMT 1.6.0, mckenzie2025 is NOT yet shipped upstream.\n' ...
         'To reproduce this block, drop mckenzie2025.m (from Thomas McKenzie)\n' ...
         'into external/amtoolbox/<unzipped-dir>/models/.\n' ...
         '未找到 mckenzie2025，已跳过 PBC-2 计算。\n' ...
         '提示：截至 AMT 1.6.0，mckenzie2025 尚未合入主线。\n' ...
         '要复现此段，请从 Thomas McKenzie 处获取 mckenzie2025.m\n' ...
         '并放入 external/amtoolbox/<解压目录>/models/。']);
end

%% Binaural frequency response (Fig. 4) / 双耳频响（图 4）
figure;
subplot(1,2,1);
freqplot_smooth(brir_orig(:,1), fs); hold on;
freqplot_smooth(brir_new (:,1), fs);
pbaspect([1.5 1 1]); ylim([0 35]); xlim([40 20000]); title('Left');
subplot(1,2,2);
freqplot_smooth(brir_orig(:,2), fs); hold on;
freqplot_smooth(brir_new (:,2), fs);
ylabel(''); yticklabels(''); pbaspect([1.5 1 1]); ylim([0 35]); xlim([40 20000]);
legend({'Original','Extrapolated'}, 'location','southwest'); title('Right');
aa = subplot(122); aa.Position(1) = 0.5;
if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_orig_bin.pdf'); end

figure;
subplot(1,2,1);
freqplot_smooth(brir_tar(:,1), fs); hold on;
freqplot_smooth(brir_new(:,1), fs);
pbaspect([1.5 1 1]); ylim([0 35]); xlim([40 20000]); title('Left');
subplot(1,2,2);
freqplot_smooth(brir_tar(:,2), fs); hold on;
freqplot_smooth(brir_new(:,2), fs);
ylabel(''); yticklabels(''); pbaspect([1.5 1 1]); ylim([0 35]); xlim([40 20000]);
legend({'Target','Extrapolated'}, 'location','southwest'); title('Right');
aa = subplot(122); aa.Position(1) = 0.5;
if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_tar_bin.pdf'); end

%% Horizontal DoA (Fig. 5) / 水平面 DoA
srir_combined = cat(3, srir_orig, srir_new, srir_tar);
figure;
[~, ~, P_pwd_h] = plot_doa_horiz(srir_combined, fs); %#ok<ASGLU>
yticklabels({'Original','Extrapolated','Target'});
set(gca, 'FontSize', 12);
pbaspect([3 1 1]);
if save_figs; exportgraphics(gcf, 'srir_extrap_hor_doa.pdf'); end

%% Spherical DoA heatmaps (Fig. 6) / 三维 DoA 热度图
figure;
[P_pwd, ~, ~, grid_dirs] = get_pwd(srir_orig, fs);
heatmap_plot(rad2deg(grid_dirs(:,1)), rad2deg(grid_dirs(:,2)), P_pwd);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
k = colorbar; xlabel(k, 'Normalised power (dB)');
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_original.pdf'); end

figure;
[P_pwd, ~, ~, grid_dirs] = get_pwd(srir_new, fs);
heatmap_plot(rad2deg(grid_dirs(:,1)), rad2deg(grid_dirs(:,2)), P_pwd);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
k = colorbar; xlabel(k, 'Normalised power (dB)');
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_extrap.pdf'); end

figure;
[P_pwd, ~, ~, grid_dirs] = get_pwd(srir_tar, fs);
heatmap_plot(rad2deg(grid_dirs(:,1)), rad2deg(grid_dirs(:,2)), P_pwd);
colormap(flipud(bone)); clim([0 3]); set(gca, 'FontSize', 11);
k = colorbar; xlabel(k, 'Normalised power (dB)');
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_target.pdf'); end

%% Audio export / 参考音频导出
if flag_render_binaural
    audioOutDir = fullfile(repoRoot, 'audio_out');
    if ~isfolder(audioOutDir); mkdir(audioOutDir); end

    speechFile = fullfile(repoRoot, 'audio_samples', 'speechdirectsound_48.wav');
    if exist(speechFile, 'file') ~= 2
        % fall back to any short dry speech sample in audio_samples / 任意 dry 样本
        cand = dir(fullfile(repoRoot, 'audio_samples', '*.wav'));
        if isempty(cand)
            warning('No dry-source audio in audio_samples/, skipping audio export.\n未在 audio_samples 找到 wav 文件，跳过音频导出。');
        else
            speechFile = fullfile(cand(1).folder, cand(1).name);
        end
    end

    if exist(speechFile, 'file') == 2
        [song_Dry, fs_src] = audioread(speechFile);
        if size(song_Dry, 2) == 2
            song_Dry = mean(song_Dry, 2);
        end
        if fs_src ~= fs
            song_Dry = resample(song_Dry, fs, fs_src);
        end

        binSong_originalIR = [conv(song_Dry, brir_orig(:,1)), conv(song_Dry, brir_orig(:,2))];
        binSong_generateIR = [conv(song_Dry, brir_new (:,1)), conv(song_Dry, brir_new (:,2))];
        binSong_recordIR   = [conv(song_Dry, brir_tar (:,1)), conv(song_Dry, brir_tar (:,2))];

        nrm = max(abs(binSong_originalIR(:)));
        if nrm == 0; nrm = 1; end

        base = sprintf('S%dL%d_to_S%dL%d', ...
            idx_LS_src_orig, idx_LS_rec_orig, idx_LS_src_tar, idx_LS_rec_tar);
        audiowrite(fullfile(audioOutDir, [base '_Original.wav']),      binSong_originalIR/nrm, fs);
        audiowrite(fullfile(audioOutDir, [base '_Extrapolation.wav']), binSong_generateIR/nrm, fs);
        audiowrite(fullfile(audioOutDir, [base '_Target.wav']),        binSong_recordIR  /nrm, fs);
        fprintf('Audio written to / 音频已导出: %s\n', audioOutDir);
    end
end

%% ===================== local functions / 局部函数 ========================

function [locD, ValD, pks, lcs] = findDirectSound(ir)
% Find the direct-sound peak after a 5 kHz high-pass.
% 用 5 kHz 高通滤波器后找到直达声峰位置与幅度。
highPassFilterFreq = 5000;
fs = 48000;
[~, filtHi, ~] = ambisonic_crossover(highPassFilterFreq, fs);
ir = filter(filtHi, 1, ir);
ir = circshift(ir, -floor(length(filtHi)/2));
absir = abs(ir);
absir = absir ./ max(absir);
[pks, lcs] = findpeaks(absir, 'MinPeakDistance', 50, 'MinPeakHeight', 0.05*max(absir));
[~, ix_max] = max(pks);
locD = lcs(ix_max);  ValD = pks(ix_max);
end

function binIR = binSound(IR_test, SH_amb_bin_dec)
% Convolve HOA signal with binaural decoder. / 将 HOA 信号与双耳解码器卷积。
ambisonic_soundscape = IR_test.';
binaural_ambisonic_render = zeros( ...
    length(SH_amb_bin_dec(1,:,1)) + length(ambisonic_soundscape) - 1, ...
    length(SH_amb_bin_dec(1,1,:)));
for i = 1:length(SH_amb_bin_dec(:,1,1))
    binaural_ambisonic_render(:,1) = binaural_ambisonic_render(:,1) + ...
        conv(SH_amb_bin_dec(i,:,1), ambisonic_soundscape(i,:))';
    binaural_ambisonic_render(:,2) = binaural_ambisonic_render(:,2) + ...
        conv(SH_amb_bin_dec(i,:,2), ambisonic_soundscape(i,:))';
end
binIR = binaural_ambisonic_render;
end

function [dist_imgsrc, delay_imgsrc, dir_imgsrc, dist_DS_imgsrc, delay_DS_imgsrc, dir_DS_imgsrc] = ...
        imgsrc_calculate(coord_rec, coord_src, ord_imgsrc, dim_room, fs)
% Calculate image sources for a shoebox room of given dimensions.
% 长方体房间内计算镜像源（坐标 / 距离 / 到达时间 / 到达方向）。
c = 343; % speed of sound (m/s)
h = figure;
plotRoom(dim_room, coord_rec, coord_src, h);
hold on;
Lx = dim_room(1); Ly = dim_room(2); Lz = dim_room(3);
x = coord_src(1); y = coord_src(2); z = coord_src(3);

if ord_imgsrc == 1
    xyz_src = [x -y -z; -x  y -z; -x -y  z; ...
               x -y -z; -x  y -z; -x -y  z].';
    xyz_room = [0 0 0; 0 0 0; 0 0 0; ...
                2*Lx 0 0; 0 2*Ly 0; 0 0 2*Lz].';
    coords_imgsrc = xyz_room - xyz_src;
    for kk = 1:size(xyz_src, 2)
        plot3(coords_imgsrc(1,kk), coords_imgsrc(2,kk), coords_imgsrc(3,kk), 'g*');
    end
else
    xyz_src = [-x -y -z; -x -y  z; -x  y -z; -x  y  z; ...
                x -y -z;  x -y  z;  x  y -z;  x  y  z].';
    nVect = -(ord_imgsrc-1):(ord_imgsrc-1);
    lVect = nVect;  mVect = nVect;
    coords_imgsrc_tot = [];
    xyz_diff_lim = [min(nVect)*2*Lx; min(lVect)*2*Ly; min(mVect)*2*Lz];
    for n = nVect
        for l = lVect
            for m = mVect
                xyz_diff = [n*2*Lx; l*2*Ly; m*2*Lz];
                coords_imgsrc = xyz_diff - xyz_src;
                delete_imgsrc = [];  i2 = 0;
                for i1 = 1:size(xyz_src, 2)
                    coord_imgsr = coords_imgsrc(:, i1);
                    if ~isempty(coord_imgsr(coord_imgsr < xyz_diff_lim))
                        i2 = i2 + 1;  delete_imgsrc(i2) = i1; %#ok<AGROW>
                    end
                end
                coords_imgsrc(:, delete_imgsrc) = [];
                coords_imgsrc_tot = [coords_imgsrc_tot coords_imgsrc]; %#ok<AGROW>
                for kk = 1:size(coords_imgsrc, 2)
                    plot3(coords_imgsrc(1,kk), coords_imgsrc(2,kk), coords_imgsrc(3,kk), 'k*');
                    plot3([coords_imgsrc(1,kk); coord_rec(1)], ...
                          [coords_imgsrc(2,kk); coord_rec(2)], ...
                          [coords_imgsrc(3,kk); coord_rec(3)], ...
                          'Color', [0.85 0.85 0.85]);
                end
            end
        end
    end
    coords_imgsrc = coords_imgsrc_tot;
end

plot3(coord_src(1), coord_src(2), coord_src(3), '^', 'LineWidth', 2, 'MarkerEdgeColor', [0.466 0.674 0.188]);
plot3(coord_rec(1), coord_rec(2), coord_rec(3), 'o', 'LineWidth', 2, 'MarkerEdgeColor', [0.85 0.325 0.098]);
view(0, 90);

dist_imgsrc  = zeros(1, size(coords_imgsrc, 2));
delay_imgsrc = dist_imgsrc;
azimuth      = dist_imgsrc;
elevation    = dist_imgsrc;
for i = 1:size(coords_imgsrc, 2)
    coord_imgsrc = coords_imgsrc(:, i);
    dist_imgsrc(i)  = norm(coord_imgsrc(:) - coord_rec(:), 2);
    delay_imgsrc(i) = round((fs/c) * dist_imgsrc(i));
    xyz_diff = coord_imgsrc - coord_rec';
    hyp = sqrt(xyz_diff(1)^2 + xyz_diff(2)^2);
    elevation(i) = atan(xyz_diff(3) / (hyp + eps));
    azimuth(i)   = atan2(xyz_diff(2), xyz_diff(1));
end
dir_imgsrc = [azimuth; elevation] * 180/pi;

% Direct sound / 直达声
dist_DS_imgsrc  = norm(coord_src(:) - coord_rec(:), 2);
delay_DS_imgsrc = round((fs/c) * dist_DS_imgsrc);
xyz_diff = coord_src' - coord_rec';
hyp = sqrt(xyz_diff(1)^2 + xyz_diff(2)^2);
elevation_DS = atan(xyz_diff(3) / (hyp + eps));
azimuth_DS   = atan2(xyz_diff(2), xyz_diff(1));
dir_DS_imgsrc = [azimuth_DS; elevation_DS] * 180/pi;
end

function [srir_arrival, srir_arrival_removed] = srir_arrival_remove( ...
        srir_in, srir_unAlt, fs, time_arrival_samp, arrival_dur_ms) %#ok<INUSL>
% Extract a windowed arrival and zero it out in the parent SRIR.
% 用窗函数在原 SRIR 上抠出一段到达分量，同时返回剔除后的余量。
arrival_window_samp = [5 10];  % linear ramps to avoid clicks / 防止"咔"声
arrival_window = ones(arrival_dur_ms/1000*fs, size(srir_in, 2));
arrival_window(1:arrival_window_samp(1), :) = ...
    repmat(linspace(0, 1, arrival_window_samp(1))', 1, size(srir_in, 2));
arrival_window(end - arrival_window_samp(2) + 1:end, :) = ...
    repmat(linspace(1, 0, arrival_window_samp(2))', 1, size(srir_in, 2));
srir_extract_window = srir_in * 0;
srir_extract_window(time_arrival_samp+1 : time_arrival_samp + arrival_dur_ms/1000*fs, :) = arrival_window;
srir_arrival_removed = srir_in .* (1 - srir_extract_window);
srir_arrival = srir_in(time_arrival_samp+1 : time_arrival_samp + arrival_dur_ms/1000*fs, :) .* arrival_window;
end

function [doa_est, normalized_doa_est_P_dB, P_pwd] = plot_doa_horiz(srir, fs)
% Plot horizontal-plane DoA over the first 0.3 s of a stack of SRIRs.
% 在多份 SRIR 的前 0.3 s 内估计水平面 DoA 并绘制散点。
res_deg_azi = 2;  res_deg_ele = 90;
order   = sqrt(size(srir, 2)) - 1;
nSrc    = 7;
numSamps = fs * 0.3;
highPassFilterFreq = 2000;  kappa = 40;
grid_dirs = grid2dirs(res_deg_azi, res_deg_ele, 0, 0);
grid_dirs = grid_dirs(2:end-1, :);  % drop poles
P_src = diag(ones(numSamps, 1));
[~, filtHi, ~] = ambisonic_crossover(highPassFilterFreq, fs);
doa_est   = zeros(nSrc, 2, size(srir, 3));
doa_est_P = zeros(nSrc, size(srir, 3));
P_pwd     = zeros(size(grid_dirs, 1), size(srir, 3));

for i = 1:size(srir, 3)
    Y_src = filter(filtHi, 1, srir(1:numSamps, :, i));
    stVec  = Y_src';
    sphCOV = stVec * P_src * stVec' + 1 * eye((order+1)^2) / (4*pi);
    [P_pwd(:, i), est_dirs_pwd, est_dirs_P] = sphPWDmap_compat(sphCOV, grid_dirs, nSrc, kappa);
    est_dirs_pwd = est_dirs_pwd * 180/pi;
    est_dirs_pwd(est_dirs_pwd < -178) = est_dirs_pwd(est_dirs_pwd < -178) + 360;
    doa_est(:, :, i) = est_dirs_pwd;
    doa_est_P(:, i) = est_dirs_P;
end

normalized_doa_est_P    = doa_est_P;
normalized_doa_est_P_dB = mag2db(normalized_doa_est_P) / 2;
plot_thresh = min(normalized_doa_est_P_dB(:));
normalized_doa_est_P_dB(normalized_doa_est_P_dB < plot_thresh) = plot_thresh;
cmap = flip(parula(256));
c    = cmap(20:end, :);
doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin', plot_thresh);
cspace = linspace(0, 1, size(c, 1));
doa_color = zeros([size(doa_01), 3]);
for i = 1:size(srir, 3)
    doa_color(:, :, i) = interp1(cspace, c(:, i), doa_01);
end
hold off;
for i = nSrc:-1:1
    scatter(squeeze(doa_est(i,1,:)), 1:size(srir, 3), 150, ...
            squeeze(doa_color(i, :, :)), 'x', 'LineWidth', 3);
    hold on;
end
xlabel('Azimuth (°)'); ylabel('');
xlim([-178 182]); xticks(-180:45:180);
colormap(c);
k = colorbar; xlabel(k, 'Normalised power (dB)');
clim([plot_thresh 0]);
set(gcf, 'Color', 'w'); pbaspect([1.7 1 1]);
box on; grid on; set(gca, 'FontSize', 16);
ylim([0.5 size(srir, 3) + 0.5]); yticks(1:size(srir, 3));
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');
end

function doa_deg = doa_iv(srir)
% Intensity-vector DoA over a (mostly) impulsive segment of an HOA signal.
% 用 W·[X,Y,Z] 强度向量法估计 HOA 信号的 DoA（输出度数）。
doa_samp = srir(:, 1) .* [srir(:, 4), srir(:, 2), srir(:, 3)];
idxNonZero = find(sum(doa_samp ~= 0, 2));
if isempty(idxNonZero)
    doa = mean(doa_samp, 1);
else
    doa = mean(doa_samp(idxNonZero, :), 1);
    doa = doa ./ vecnorm(doa, 2, 2);
end
[doa1(1,1), doa1(2,1), ~] = cart2sph(doa(1), doa(2), doa(3));
doa_deg = rad2deg(doa1);
end

function [P_pwd, doa_est, doa_est_P, grid_dirs] = get_pwd(srir, fs)
% Plane-wave decomposition over a 2°×2° sphere grid for a single SRIR.
% 对单条 SRIR 做 2°×2° 球面网格上的平面波分解。
degreeResolution = 2;
order   = sqrt(size(srir, 2)) - 1;
nSrc    = 7;
numSamps = fs * 0.3;
highPassFilterFreq = 3000; kappa = 40;
grid_dirs = grid2dirs(degreeResolution, degreeResolution, 0, 0);
P_src = diag(ones(numSamps, 1));
[~, filtHi, ~] = ambisonic_crossover(highPassFilterFreq, fs);
doa_est = zeros(nSrc, 2);
doa_est_P = zeros(nSrc, 1);
P_pwd = zeros(size(grid_dirs, 1), 1);
Y_src = filter(filtHi, 1, srir(1:numSamps, :));
stVec  = Y_src';
sphCOV = stVec * P_src * stVec' + 1 * eye((order+1)^2) / (4*pi);
[P_pwd(:), est_dirs_pwd, est_dirs_P] = sphPWDmap_compat(sphCOV, grid_dirs, nSrc, kappa);
est_dirs_pwd = est_dirs_pwd * 180/pi;
est_dirs_pwd(est_dirs_pwd < -170) = est_dirs_pwd(est_dirs_pwd < -170) + 360;
doa_est(:, :) = est_dirs_pwd;
doa_est_P(:) = est_dirs_P;
end

function plotRoom(roomDimensions, receiverCoord, sourceCoord, figHandle)
% Draw the shoebox room with source and receiver markers.
% 绘制鞋盒房间并标出声源、接收点。基于 MATLAB 官方示例改编。
figure(figHandle);
X = [0; roomDimensions(1); roomDimensions(1); 0; 0];
Y = [0; 0; roomDimensions(2); roomDimensions(2); 0];
Z = [0; 0; 0; 0; 0];
figure; hold on;
plot3(X, Y, Z, 'k', 'LineWidth', 1.5);
plot3(X, Y, Z + roomDimensions(3), 'k', 'LineWidth', 1.5);
set(gca, 'View', [-28, 35]);
for k = 1:length(X)-1
    plot3([X(k); X(k)], [Y(k); Y(k)], [0; roomDimensions(3)], 'k', 'LineWidth', 1.5);
end
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
plot3(sourceCoord(1),   sourceCoord(2),   sourceCoord(3),   'bx', 'LineWidth', 2);
plot3(receiverCoord(1), receiverCoord(2), receiverCoord(3), 'ro', 'LineWidth', 2);
end
