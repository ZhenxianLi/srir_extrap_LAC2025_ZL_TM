% Subspace Decomposition / 子空间分解
% Cut Direct Sound / 剪切直达声
% Scale Direct Sound / 缩放直达声
% Rotate Direct Sound / 旋转直达声
% Rejoin / 重新拼接

% From Single Recording / 单次录音数据处理
% --------------------------------------------------------------------------------------------------------------------------
%% Load and Plot SOFA File / 加载并绘制SOFA文件

% Includes scaling early reflections and simulating spatial audio responses 
% at different room positions by modifying the direct sound arrival time and amplitude.
% 包括对早期反射进行缩放并通过改变直达声到达时间与幅度来模拟在房间内不同位置处的空间音频响应。

clc;
clear all;
close all;

% Start SOFA and Load Function Library / 启动SOFA并加载函数库
% Loads SOFA API and other required libraries for handling Higher-Order Ambisonics and SOFA format impulse responses.
% 加载SOFA相关API以及其他所需函数库，用于处理高阶Ambisonics和SOFA格式的脉冲响应

% Get Parent Folder Path / 获取父文件夹路径
% Retrieves the parent folder of the current working directory for adding dependency paths.
parentFolder = fileparts(pwd);

addpath(fullfile(parentFolder, 'API_MO'));
addpath(fullfile(parentFolder, 'Spherical-Harmonic-Transform-master'));
addpath(fullfile(parentFolder, 'Higher-Order-Ambisonics-master'));
addpath(genpath(fullfile(parentFolder, 'binaural_ambisonic_preprocessing-main')));
addpath(fullfile(parentFolder, 'audio_samples'));
addpath(fullfile(parentFolder, 'SRIR-Subspace-Decomposition-master'));
SOFAstart;

% Run to Load Decoder / 加载解码器
% If Ambisonics configuration is not loaded, run load_ambisonic_configuration to obtain the Ambisonics decoder.
if ~exist('aio_flag')
    load_ambisonic_configuration
end

% Load Variable Acoustics 6DoF Dataset (Most Reverberant Set) / 加载可变声学6自由度数据集（最混响版本）
% Loads the 6DoF SRIRs data stored in the '6dof_SRIRs_eigenmike_SH' folder.
parentFolder = fileparts(pwd);
irPath1 = fullfile(parentFolder, '6dof_SRIRs_eigenmike_SH/');
irName1 = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';

% Download the SRIR.sofa from Github Release / 从Github下载SRIR.sofa文件
% If the SOFA file is not available locally, download it from Github.
url = 'https://github.com/ZhenxianLi/ZHENXIAN_SRIRextrapolation/releases/download/6dof_SRIRs_eigenmike_SH/6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';
% Construct the full path for the saved file / 构造保存文件的完整路径
outputFile = fullfile(irPath1, irName1);
% Check if the file already exists to avoid duplicate downloads / 检查文件是否存在，避免重复下载
if ~exist(outputFile, 'file')
    fprintf('start Download the file: %s\n', irName1);
    try
        % download to pass
        websave(outputFile, url);
        fprintf('downloadSucess: %s\n', outputFile);
    catch ME
        fprintf('downloadFail: %s\n', ME.message);
        return; % Terminate execution / 终止执行
    end
else
    fprintf('File already exists, skipping download: %s\n', outputFile);
end

% Load the SOFA file using SOFA API / 使用SOFAAPI加载SOFA文件
sofa1 = SOFAload([irPath1,irName1]);
fs = sofa1.Data.SamplingRate;
% SOFAplotGeometry(sofa1);   % Plot the source and listener positions / 绘制声源与监听点位置
close all;

%% Check Information / 检查信息
% Display metadata and dimensions from the SOFA file.
% 显示SOFA文件中的元数据和尺寸信息。
sofa1.GLOBAL_RoomType;
sofa1.GLOBAL_SOFAConventionsVersion;
sofa1.ListenerPosition;
sofa1.SourcePosition;
sofa1.EmitterPosition;
size(sofa1.Data);

% Check the first reflection index manually / 手动检查第一早期反射索引
% In each IR group, the manually specified index for early reflection arrival, which can be obtained via peak detection.
firstReflectionIndex = [487;697;674;528;579;621;606;515;616;726;499;658;861;543;742;939;971;1146;1071;900;1171];

%% Calculate Sound Speed / 计算声速
% Get IR dimensions and select only N=15 (only keep 15 front-facing directions).
% IROrder: Number of channels corresponding to Higher-Order Ambisonics order.
% L: Length of each IR.
sizeIR = size(sofa1.Data.IR);
N = sizeIR(1);
N = 15;     % only facing=15 / 只取扬声器正面朝向麦克风的方向的组合，即前15个
IROrder = sizeIR(2); % How many orders contained in the SRIR recording / 高阶Ambisonics通道数量
L = sizeIR(3);        % Sample length of the signal / 单条IR的采样长度

% Find direct sound and first early reflection arrival time / 寻找直达声和第一早期反射到达时间
% Estimate propagation speed based on the direct sound arrival index and corresponding distance.
speedSum = 0;
distance = zeros(N,1);
directSoundTime = zeros(N,1);
directSoundValure = zeros(N,1);
for n = 1:N
    IR_n = squeeze(sofa1.Data.IR(n,1,:)); % n, 1st order, full length
    [directSoundIndex, directSoundValure(n), pks, lcs] = findDirectSound(IR_n);
    distance(n) = sqrt(sum((sofa1.ListenerPosition(n,:) - sofa1.SourcePosition(n,:)).^2));
    directSoundTime(n) = directSoundIndex;  % Time of direct sound (in samples) / 直达声到达时间（采样点）
    speedSum = speedSum + distance(n) / directSoundTime(n);
    speed = speedSum / n; % in samples
end

speedinSec = speed * fs;

% Plot all N sounds / 绘制所有N个声
% If pic1==1, visualize the direct sound positions and early reflection peaks.
pic1 = 0;
if pic1
    figure;
    for n = 1:N
        IR_n = squeeze(sofa1.Data.IR(n,1,:)); % n, 1st order, full length
        subplot(7,3,n); % 3 speakers; 7 mics / 3个扬声器；7个麦克风
        plot(abs(IR_n));
        hold on;
        xline(directSoundIndex, '--g');
        xline(firstReflectionIndex(n), '--g');
        ylim([0,0.5]);
        xlim([0,4000]);
        hold on;
        plot(lcs, pks);
        title(num2str(n));
        hold on;
    end
    hold off;
    % end for loop
end

% Generate New IR / 生成新的IR
%% Subspace Decomposition and Cut / 子空间分解与剪切
% Begin generating the target SRIR. Use the 9th SRIR (S3L3, room center) as the reference,
% and scale according to the distance and direct sound time delay/amplitude differences.
% 以第9条SRIR（S3L3，房间中心）作为参考，通过计算目标位置与声源之间的距离，以及直达声时延和振幅的缩放来生成目标SRIR。

% Origin Position / 原始位置
% Speaker: originSpeaker, Mic: originListener (used as reference)
originSpeaker = 1; % (1-3)
originListener = 1; % (1-7)
% Target Position / 目标位置
% Speaker: targetSpeaker, Mic: targetListener (combination to be generated)
targetSpeaker = 1; % (1-3)
targetListener = 5; % (1-7)

% Calculate record indices for origin and target positions / 计算原始和目标位置的录音索引
originRecordIndex = originListener * 3 - 3 + originSpeaker; % Index for origin SRIR
targetRecordIndex = targetListener * 3 - 3 + targetSpeaker;  % Index for target SRIR

% Set cut parameters / 设置剪切参数
directSoundCutLeft = 50;
directSoundCutRight = 300;
earlyRefCutLength = 400;

% Read IR and positions / 读取IR和位置信息
% IR_original is the reference IR (9th recording) / IR_original 为参考IR（第9条）
IR_original = sofa1.Data.IR(originRecordIndex,:,:);

SourcerPoint_Origin = sofa1.SourcePosition(originRecordIndex,:);
ListenerPoint_Origin = sofa1.ListenerPosition(originRecordIndex,:);

SourcerPoint_generate = sofa1.SourcePosition(targetRecordIndex,:);
ListenerPoint_generate = sofa1.ListenerPosition(targetRecordIndex,:); % Recording for target position / 目标位置的录音
arrival_time_original = directSoundTime(originRecordIndex); % in samples / 采样点

%% Decomposition of Direct and Residual / 分解直达声与剩余部分

srir_for_decomp = squeeze(IR_original).';
fs = fs; % Sampling frequency remains unchanged / 采样频率不变

% Parameters for subspace decomposition / 子空间分解参数
% See function header of srirSubspaceDecomp for details.
blockLenSmp = 32;
hopSizeSmp = blockLenSmp / 8;
kappa = 3;
numBlocksGsvSumAvg = 32;
residualEstimateLengthMs = 20;
decompositionTimeLimitMs = 100;
numBlocksSmoothThresh = 1;

[dirSrir, resSrir, numDirSubspaceComponents, gsvs, detectionThreshold, gsvSum, avgGsvSum] = ...
            srirSubspaceDecomp(srir_for_decomp, fs, blockLenSmp, hopSizeSmp, kappa, numBlocksGsvSumAvg, residualEstimateLengthMs, ...
                               decompositionTimeLimitMs, numBlocksSmoothThresh);

%% Plot the Decomposition Result / 绘制分解结果
plot_decomposition_result = 1;
if plot_decomposition_result
    srirLen = size(srir_for_decomp,1);
    t = linspace(0, srirLen/fs - 1/fs, srirLen).';
    tBlocks = (0:hopSizeSmp:size(gsvs,1)*hopSizeSmp - hopSizeSmp) / fs;
    
    % Plot cumulative sum of direct and residual parts
    fig_decomposition_sum = figure('Name', 'Decomposition Sum - All Channels');
    hold on;
    plot(t*1000, sum(abs(dirSrir), 2), 'LineWidth', 2);
    plot(t*1000, sum(abs(resSrir), 2), 'LineWidth', 2);
    xlim([0,100]);
    xlabel('$t$ (ms)', 'Interpreter', 'latex');
    ylabel('$\| \cdot \|$ (dB)', 'Interpreter', 'latex');
    legend({'$\mathbf{x}_\mathrm{d}(t)$', '$\mathbf{x}_\mathrm{r}(t)$'}, 'Interpreter', 'latex');
    grid on;
    
    % Plot GSV (Generalized Singular Values) cumulative sums
    numChannels = size(srir_for_decomp,2);
    cumsumGSVsColors = copper(numChannels);
    cumsumGSVs = cumsum(gsvs,2);
    
    fig_decomposition2 = figure('Name', 'Decomposition GSV');
    hold on;
    for ii = 1:numChannels
        if ii == numChannels
            hGSVCumsum = plot(tBlocks*1000, cumsumGSVs(:,ii), 'Color', cumsumGSVsColors(ii,:), 'LineWidth', 1.5);
        else
            plot(tBlocks*1000, cumsumGSVs(:,ii)*numChannels/ii, 'Color', cumsumGSVsColors(ii,:), 'LineWidth', 1.5);
        end
    end
    hAvgGSVs = plot(tBlocks*1000, avgGsvSum, 'k:', 'LineWidth', 1.5);
    hDetectionThresh = plot(tBlocks*1000, detectionThreshold, 'k', 'LineWidth', 1.5);
    grid on;
    xlabel('$t$ (ms)', 'Interpreter', 'latex');
    xlim([0,100]);
    ylim([0,10]);
    legend([hGSVCumsum, hAvgGSVs, hDetectionThresh], ...
        {'Cumulative Sums of GSVs', 'Subspace Component Threshold', 'Detection Threshold'}, ...
        'Interpreter', 'latex');
end

% Transpose matrices to maintain original code orientation / 倒置矩阵以保持原始代码形式一致
dirSrir = dirSrir.';
resSrir = resSrir.';

%% Scale the Direct Part / 缩放直达声部分
% Fit distance-amplitude model using RMS values.
% For Nface=15, compute the RMS of the direct sound across 15 front-facing directions
% to fit an exponential decay model based on distance.
Nface = 15; % Only facing 15 directions / 只取15个正面方向
rmsDirectSoundValure = zeros(Nface,1);

for n = 1:Nface
    readIR_forFit = sofa1.Data.IR(n,:,:);
    arrival_time_forFit = firstReflectionIndex(n);
    ds = readIR_forFit(:,1,arrival_time_forFit - directSoundCutLeft : arrival_time_forFit + directSoundCutRight);
    rmsDirectSoundValure(n) = rms(squeeze(ds));
end

% Define a function handle for fitting the exponential decay model between distance and RMS
% 用于拟合距离与RMS之间的指数衰减模型
fun = @(x) sseval(x, distance(1:Nface,1), rmsDirectSoundValure);
x0 = rand(2,1);
bestx = fminsearch(fun, x0);

% Visualize the fitting result / 可视化拟合结果
pic2 = 1;
if pic2
    fig_distance_scale_fit = figure('Name', 'Distance Scale Fit');
    scatter(distance, rmsDirectSoundValure);
    for n = 1:N
        text(distance(n), rmsDirectSoundValure(n), num2str(n));
    end
    hold on;
    dislist = 0:0.2:10;
    plot(dislist, bestx(1)*exp(-bestx(2)*dislist));
    xlabel('Distance (m)');
    ylabel('RMS Amplitude of Direct Sound');
end

% Compute time delay and amplitude for target position / 计算目标位置的直达声时延和幅度
distance_generate = sqrt(sum((SourcerPoint_generate - ListenerPoint_generate).^2));
distance_origianl = sqrt(sum((SourcerPoint_Origin - ListenerPoint_Origin).^2));
arrival_time_generate = floor(distance_generate / speed); % in samples / 采样点
Am_generate = bestx(1) * exp(-bestx(2) * distance_generate);
Am_original = bestx(1) * exp(-bestx(2) * distance_origianl);

% Extract the direct sound segment from dirSrir / 从dirSrir中提取直达声部分
directSound = dirSrir(:, arrival_time_original - directSoundCutLeft : arrival_time_original + directSoundCutRight);
% Scale the direct sound segment according to amplitude ratio / 根据幅度比例缩放直达声部分
directSound = directSound * (Am_generate / Am_original); % scale

%% Image Source Method to Get Early Reflection / 使用镜像源法获取早期反射
% Example: Define room dimensions, source and microphone positions, with speed = speedinSec m/s.
roomSize = [7.87, 5.75, 2.91];       % Room dimensions (L, W, H) / 房间尺寸 (长, 宽, 高)
sourcePos = sofa1.SourcePosition(originRecordIndex,:);  % Source position / 声源坐标
micPos    = sofa1.ListenerPosition(originRecordIndex,:); % Microphone position / 麦克风坐标
maxRef    = 2;             % Maximum number of reflections to compute / 计算的最大镜像次数

% Call the function / 调用函数
T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxRef, speedinSec, sofa1); % Image source method and plot
T.ArrTime

figure(fig_decomposition_sum);
hold on;
for k = 1:length(T.ArrTime)
    sumRefOrder = abs(T.Nx(k)) + abs(T.Ny(k)) + abs(T.Nz(k));
    if sumRefOrder == 0
        xline(T.ArrTime(k)*1000, 'r--', ['Nxyz= ' num2str(T.Nx(k)) ' ' num2str(T.Ny(k)) ' ' num2str(T.Nz(k))], 'DisplayName', 'Direct Sound');
    elseif sumRefOrder == 1
        xline(T.ArrTime(k)*1000, 'g--', ['Nxyz= ' num2str(T.Nx(k)) ' ' num2str(T.Ny(k)) ' ' num2str(T.Nz(k))], 'DisplayName', 'First Reflection');
    elseif sumRefOrder >= 1
        xline(T.ArrTime(k)*1000, 'b--', ['Nxyz= ' num2str(T.Nx(k)) ' ' num2str(T.Ny(k)) ' ' num2str(T.Nz(k))], 'DisplayName', 'Second Reflection');
    end
    % Alternative commented lines omitted.
end
hold off;

%% Re-join / 重新拼接
% A: Section before direct sound (padded with zeros)
% B: Direct sound slice
% C: Early reflections
% D: Late reverberation
% Note: earlyRefCutLength=400 ensures inclusion of early reflections.
A = zeros(25, arrival_time_generate - directSoundCutLeft - 1);
B = squeeze(directSound);
C = squeeze(IR_original(:,:, arrival_time_generate + directSoundCutRight + 1 : arrival_time_generate + directSoundCutRight + earlyRefCutLength)); % Early reflections
% C = C*(Am_generate/Am_original)/0.25; % Scale C (this line is commented out) / （此行缩放被注释掉）
D = squeeze(IR_original(:,:, arrival_time_generate + directSoundCutRight + earlyRefCutLength + 1 : end)); % Reverberation

IR_scale = [A, B, C, D];

%% Plot and Listen / 绘制并试听
% Optionally compare the processed signal with the original signal.
pic3 = 0;
if pic3
    % sound(IR_scale(1,:), fs);
    figure;
    
    subplot(3,1,1);
    plot(abs(squeeze(IR_original(:,1,:))));
    legend('IR Original');
    xlim([0 3000]);
    ylim([0 0.5]);
    
    subplot(3,1,2);
    plot(abs(squeeze(IR_scale(1,:))));
    legend('IR Generated');
    xlim([0 3000]);
    ylim([0 0.5]);
    
    subplot(3,1,3);
    plot(abs(squeeze(sofa1.Data.IR(targetRecordIndex,1,:))));
    legend('IR Recorded');
    xlim([0 3000]);
    ylim([0 0.5]);
end

%% Rotate / 旋转
% Rotate the HOA signal by computing the deviation in spherical coordinates (azimuth, elevation)
% between (SourcerPoint_Origin - ListenerPoint_Origin) and (SourcerPoint_generate - ListenerPoint_generate).
% 通过计算两个位置在球坐标系下（方位角、仰角）的偏差来实现HOA信号旋转
% (SourcerPoint_Origin - ListenerPoint_Origin) 与 (SourcerPoint_generate - ListenerPoint_generate)
 
% Compute spherical coordinates for both positions / 计算两个位置的球坐标
xyz1 = SourcerPoint_Origin - ListenerPoint_Origin;
xyz2 = SourcerPoint_generate - ListenerPoint_generate;
[az1, el1, r1] = cart2sph(xyz1(1), xyz1(2), xyz1(3));
[az2, el2, r2] = cart2sph(xyz2(1), xyz2(2), xyz2(3));
yaw = (az2 - az1) / pi * 180;
pitch = (el2 - el1) / pi * 180;
roll = 0;

% Reshape HOA signal / 重塑HOA信号
% Transpose B (direct sound part) to match the input format of rotateHOA_N3D function.
hoasig = B.'; % Reshape direct sound

% Rotate / 旋转
% Use rotateHOA_N3D to perform HOA signal rotation.
hoasig_rot = rotateHOA_N3D(hoasig, yaw, pitch, roll);
B_rot = hoasig_rot.';

% Re-join / 重新拼接
% Combine the rotated direct sound part with early reflections and late reverberation.
IR_rot = [A, B_rot, C, D];

%% Plot and Listen (Post-Rotation) / 绘制并试听（旋转后）
% Compare amplitude differences between original, rotated, and actual recorded IR.
pic3 = 1;
if pic3
    % sound(IR_scale(1,:), fs);
    figure;
    
    subplot(3,1,1);
    plot(abs(squeeze(IR_original(:,1,:))));
    legend(['IR Original ' num2str(3) num2str(3)]);
    xlim([0 8000]);
    ylim([0 0.5]);
    
    subplot(3,1,2);
    plot(abs(squeeze(IR_rot(1,:))));
    legend('IR Rotated');
    xlim([0 8000]);
    ylim([0 0.5]);
    
    subplot(3,1,3);
    plot(abs(squeeze(sofa1.Data.IR(targetRecordIndex,1,:))));
    legend(['IR Recorded ' num2str(targetSpeaker) num2str(targetListener)]);
    xlim([0 8000]);
    ylim([0 0.5]);
end

%% Render Binaural and Play / 双耳渲染并播放
% Convolve HOA signal with binaural decoder to obtain binaural signal for listening.
 
% Play all records / 播放所有录音
% If playAllRecord==1, loop through and play all recordings.
playAllRecord = 0;
if playAllRecord
    for i = 1:N/3
        speakerNum = 3;
        recordNum = i * 3 - 3 + speakerNum;
        IR_test = squeeze(sofa1.Data.IR(recordNum,:,:));
        waitingTime = zeros((i-1) * 124223, 2);
        binSound(IR_test, SH_ambisonic_binaural_decoder, Fs, waitingTime);
    end
end

% Render and play: IR_original (9), IR_generated, and IR_recorded / 渲染并播放：原始IR（第9条）、生成IR、实测IR
% Compare: 1. S3L3, 2. Generated, 3. Target recording.
doBinRenderIR = 1;
if doBinRenderIR
    IR_Record = squeeze(sofa1.Data.IR(targetRecordIndex,:,:));
    IR_original = squeeze(sofa1.Data.IR(9,:,:));
    
    binIR_original = binSound(IR_original, SH_ambisonic_binaural_decoder);
    binIR_generate = binSound(IR_rot, SH_ambisonic_binaural_decoder);
    binIR_record = binSound(IR_Record, SH_ambisonic_binaural_decoder);
    gaptime = zeros(0.3 * fs, 2);
    soundsc([binIR_original; gaptime; binIR_generate; gaptime; binIR_record], Fs);
    % soundsc([binIR_generate; gaptime; binIR_record], Fs); % Only last two sounds / 仅播放后两个信号
end

% Play convolved dry signal with IRs / 用干声信号与IR卷积播放
% Compare: 1. S3L3, 2. Generated, 3. Recorded.
doBinRenderSONG = 0;
if doBinRenderSONG
    % song_Dry = audioread('speechdirectsound_48.wav');
    [~, ch] = size(song_Dry);
    if ch == 2
        song_Dry = (song_Dry(:,1) + song_Dry(:,2)) / 2;
    end
    
    clear binSong_recordIR;
    clear binSong_originalIR;
    clear binSong_generateIR;
    
    for i = 1:2 % Left and Right channels / 左右声道
        binSong_originalIR(:,i) = conv(song_Dry, binIR_original(:,i));
        binSong_generateIR(:,i) = conv(song_Dry, binIR_generate(:,i));
        binSong_recordIR(:,i)   = conv(song_Dry, binIR_record(:,i));
    end
    
    gaptime = zeros(0.5 * fs, 2); % Gap of 0.5 seconds / 间隔0.5秒
    soundsc([binSong_originalIR; gaptime; binSong_generateIR; gaptime; binSong_recordIR], fs);
    
    if ~exist('Max')
        Max = 1;
    end
    
    audiowrite(['S' num2str(targetSpeaker) 'L' num2str(targetListener) '_Original33.wav'], binSong_originalIR/Max, fs);
    audiowrite(['S' num2str(targetSpeaker) 'L' num2str(targetListener) '_Rotate.wav'], binSong_generateIR/Max, fs);
    % audiowrite(['S' num2str(useSpeaker) 'L' num2str(useListener) '_Record.wav'], binSong_recordIR/Max, fs);
end

%% Declare Function: findDirectSound / 查找直达声位置及幅值
% This function locates the direct sound's index and amplitude in the IR
% using findpeaks with a threshold.
% 该函数利用findpeaks在设定阈值后找出直达声在IR中的索引和幅值。
function [locD, ValD, pks, lcs] = findDirectSound(ir)
    absir = abs(ir);
    noiseValuse = max(absir);
    [pks, lcs] = findpeaks(absir, "MinPeakDistance", 10, MinPeakHeight = 0.2 * noiseValuse);  % Denoising, find peak
    % Return the first peak as direct sound / 返回第一个峰值作为直达声
    locD = lcs(1);
    ValD = pks(1);
end

%% Declare Function: sseval / 计算距离-幅度拟合误差
% This function computes the sum of squared errors for the exponential decay model,
% used with fminsearch to find optimal parameters.
% 该函数计算指数衰减模型的误差和，用于fminsearch拟合最优参数A0和alpha。
function sse = sseval(x, tdata, ydata)
    A = x(1);
    lambda = x(2);
    sse = sum((ydata - A * exp(-lambda * tdata)).^2);
end

%% Declare Function: binSound / 双耳渲染函数
% Converts Higher-Order Ambisonics IR to a binaural signal by convolving with a binaural decoder.
function binIR = binSound(IR_test, SH_ambisonic_binaural_decoder)
    % Render to binaural audio and play it.
    ambisonic_soundscape = IR_test.';
    
    binaural_ambisonic_render = zeros(length(SH_ambisonic_binaural_decoder(1,:,1)) + length(ambisonic_soundscape) - 1, length(SH_ambisonic_binaural_decoder(1,1,:)));
    
    % Convolve each channel of the encoded signal with the corresponding decoder channel and sum the result.
    for i = 1:length(SH_ambisonic_binaural_decoder(:,1,1))
        binaural_ambisonic_render(:,1) = binaural_ambisonic_render(:,1) + conv(SH_ambisonic_binaural_decoder(i,:,1), ambisonic_soundscape(:,i));
        binaural_ambisonic_render(:,2) = binaural_ambisonic_render(:,2) + conv(SH_ambisonic_binaural_decoder(i,:,2), ambisonic_soundscape(:,i));
    end
    
    % Return the binaural render.
    binIR = binaural_ambisonic_render;
end

%% Declare Function: compute_DOA_ISM_with_plot / 使用镜像源法计算并可视化DOA、到达时间等信息
% T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxReflections, c)
%
% Input:
%   roomSize       : [Lx, Ly, Lz] Room dimensions in meters.
%   sourcePos      : [xs, ys, zs] Source position in meters.
%   micPos         : [xm, ym, zm] Microphone position in meters.
%   maxReflections : Maximum number of reflections (recommended 1~2).
%   c              : Speed of sound (default 343 m/s).
%
% Output:
%   T : Table containing columns:
%       Nx, Ny, Nz : Reflection counts (for x, y, z directions).
%       Ximg, Yimg, Zimg : Image source positions.
%       Distance : Distance from image source to microphone (m).
%       ArrTime  : Arrival time (s).
%       Theta    : Azimuth (radians).
%       Phi      : Elevation (radians).
%
% Example:
%   roomSize = [7.87, 5.75, 2.91];
%   sourcePos = [2, 1, 1.2];
%   micPos    = [3, 2, 1];
%   maxRef    = 2;
%   c         = 343;
%   T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxRef, c);
function T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxReflections, c, sofa1)
    % Parse room dimensions / 解析房间尺寸
    Lx = roomSize(1);
    Ly = roomSize(2);
    Lz = roomSize(3);
    
    xs = sourcePos(1); ys = sourcePos(2); zs = sourcePos(3);
    xm = micPos(1);   ym = micPos(2);   zm = micPos(3);
    
    % Initialize result storage / 初始化结果存储
    Nx_arr = [];
    Ny_arr = [];
    Nz_arr = [];
    
    Ximg_arr = [];
    Yimg_arr = [];
    Zimg_arr = [];
    
    Dist_arr = [];
    ArrTime_arr = [];
    
    Theta_arr = [];
    Phi_arr   = [];
    
    % Counters: direct, first, and second reflections / 计数器：直达声、一/二次反射
    directCount  = 0;
    singleCount  = 0;
    doubleCount  = 0;
    
    figure('Name', 'ISM Visualization', 'NumberTitle', 'off');
    hold on;
    grid on;
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Image Source Method (ISM) Visualization');
    
    % Draw room frame (optional) / 绘制房间边界（可选）
    drawRoomFrame(0, 0, 0, Lx, Ly, Lz, sofa1);
    
    % Plot microphone position / 绘制麦克风位置
    plot3(xm, ym, zm, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    text(xm, ym, zm, '  Mic', 'Color', 'k', 'FontWeight', 'bold');
    
    % Plot original source / 绘制原始声源
    plot3(xs, ys, zs, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(xs, ys, zs, '  Source', 'Color', 'r', 'FontWeight', 'bold');
    
    % Loop over reflection orders / 遍历不同镜像次数
    for nx = -maxReflections:maxReflections
        for ny = -maxReflections:maxReflections
            for nz = -maxReflections:maxReflections
                
                sumOrder = abs(nx) + abs(ny) + abs(nz);
                
                % Allow sumOrder = 0 (direct sound), 1 (first reflection), 2 (second reflection), etc.
                % 只允许 sumOrder <= maxReflections
                if sumOrder <= maxReflections
                    % Calculate image source position considering parity / 计算镜像源位置（考虑奇偶性）
                    
                    if nx >= 0
                        % For positive reflection: number of reflections on right wall = floor((nx+1)/2)
                        x_img = (-1)^nx * xs + 2 * Lx * floor((nx+1)/2);
                    else
                        % For negative reflection: number on left wall = floor(-nx/2)
                        x_img = (-1)^nx * xs - 2 * Lx * floor(-nx/2);
                    end
                    if ny >= 0
                        y_img = (-1)^ny * ys + 2 * Ly * floor((ny+1)/2);
                    else
                        y_img = (-1)^ny * ys - 2 * Ly * floor(-ny/2);
                    end
                    
                    if nz >= 0
                        z_img = (-1)^nz * zs + 2 * Lz * floor((nz+1)/2);
                    else
                        z_img = (-1)^nz * zs - 2 * Lz * floor(-nz/2);
                    end
                    
                    % Direction vector: from image source to microphone / 方向向量：镜像源到麦克风
                    d = [x_img - xm, y_img - ym, z_img - zm];
                    
                    % Calculate distance / 计算距离
                    dist = norm(d);
                    % Calculate arrival time / 计算到达时间
                    arrTime = dist / c;
                    
                    % Calculate DOA (azimuth & elevation) / 计算方位角和仰角
                    theta = atan2(d(2), d(1));   % Azimuth / 方位角
                    phi   = atan2(d(3), sqrt(d(1)^2 + d(2)^2)); % Elevation / 仰角
                    
                    % Store results / 存储结果
                    Nx_arr = [Nx_arr; nx];
                    Ny_arr = [Ny_arr; ny];
                    Nz_arr = [Nz_arr; nz];
                    
                    Ximg_arr = [Ximg_arr; x_img];
                    Yimg_arr = [Yimg_arr; y_img];
                    Zimg_arr = [Zimg_arr; z_img];
                    
                    Dist_arr = [Dist_arr; dist];
                    ArrTime_arr = [ArrTime_arr; arrTime];
                    
                    Theta_arr = [Theta_arr; theta];
                    Phi_arr   = [Phi_arr; phi];
                    
                    % Visualize based on reflection order / 根据反射次数可视化
                    if sumOrder == 0
                        directCount = directCount + 1;
                        % Mark direct sound with magenta (different from SourcePos)
                        plot3(x_img, y_img, z_img, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                        line([xm x_img], [ym y_img], [zm z_img], 'Color', 'm', 'LineStyle', '-');
                    elseif sumOrder == 1
                        singleCount = singleCount + 1;
                        % First reflection -> green
                        plot3(x_img, y_img, z_img, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
                        line([xm x_img], [ym y_img], [zm z_img], 'Color', [0.3 0.8 0.3], 'LineStyle', '--');
                    elseif sumOrder == 2
                        doubleCount = doubleCount + 1;
                        % Second reflection -> blue
                        plot3(x_img, y_img, z_img, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
                        line([xm x_img], [ym y_img], [zm z_img], 'Color', [0.3 0.3 0.8], 'LineStyle', '--');
                    else
                        % For three or more reflections, additional else-if blocks can be added.
                    end
                end
            end
        end
    end
    
    % Output the table / 输出结果表格
    T = table(Nx_arr, Ny_arr, Nz_arr, ...
        Ximg_arr, Yimg_arr, Zimg_arr, ...
        Dist_arr, ArrTime_arr, Theta_arr, Phi_arr, ...
        'VariableNames', {'Nx','Ny','Nz','Ximg','Yimg','Zimg','Distance','ArrTime','Theta','Phi'});
    
    % Display table in command window / 在命令行中显示表格
    disp(T);
    
    % Display statistics / 显示统计结果
    fprintf('\n--- Reflection Stats ---\n');
    fprintf('  Direct sound (sumOrder=0) count: %d\n', directCount);
    fprintf('  First reflection (sumOrder=1) count: %d\n', singleCount);
    fprintf('  Second reflection (sumOrder=2) count: %d\n', doubleCount);
    
    % Adjust view angle / 调整视角
    view(3);
    hold off;
end

%% Auxiliary Function: drawRoomFrame / 绘制房间外框和所有测量点
% Draws a wireframe of the room and plots all measurement points.
% 在当前坐标轴中绘制房间线框及所有测量点
function drawRoomFrame(x0, y0, z0, Lx, Ly, Lz, sofa1)
    % Define room corner points / 定义房间角点
    X = [x0, x0+Lx, x0+Lx, x0,     x0,     x0+Lx, x0+Lx, x0];
    Y = [y0, y0,     y0+Ly, y0+Ly, y0,     y0,     y0+Ly, y0+Ly];
    Z = [z0, z0,     z0,     z0,     z0+Lz, z0+Lz, z0+Lz, z0+Lz];
    edges = [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4];
    for i = 1:length(edges)-1
        plot3([X(edges(i)) X(edges(i+1))], [Y(edges(i)) Y(edges(i+1))], [Z(edges(i)) Z(edges(i+1))], 'k-', 'LineWidth', 1.0);
        hold on;
    end
    
    % Plot all measurement points / 绘制所有测量点
    for n = 1:21
        plot3(sofa1.SourcePosition(n,1), sofa1.SourcePosition(n,2), sofa1.SourcePosition(n,3), 'ko');
        hold on;
        plot3(sofa1.ListenerPosition(n,1), sofa1.ListenerPosition(n,2), sofa1.ListenerPosition(n,3), 'ko');
        hold on;
    end
end