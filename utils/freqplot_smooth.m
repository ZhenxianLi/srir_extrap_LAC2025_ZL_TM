function freqplot_smooth(x, fs, fracOct)
% FREQPLOT_SMOOTH  Fractional-octave smoothed magnitude spectrum on a log axis.
%                  对单通道时域信号做分数倍频程平滑后的幅度谱绘图。
%
% Usage / 用法:
%   freqplot_smooth(x, fs)            % default 1/3-octave smoothing / 默认 1/3 倍频程平滑
%   freqplot_smooth(x, fs, fracOct)   % e.g. 6 for 1/6-octave / 例如 6 即 1/6 倍频程
%
% Inputs / 输入:
%   x        - time-domain signal column / 单声道时间序列 (N x 1)
%   fs       - sampling rate (Hz) / 采样率
%   fracOct  - smoothing denominator (1/fracOct octave) / 平滑分母
%
% This function recreates the `freqplot_smooth` helper that the original LAC
% scripts call but did not ship. It uses a simple geometric-mean smoothing
% across fractional-octave bands and plots magnitude in dB versus log Hz.
% 该函数复刻 LAC 主脚本里调用却没有附带的 freqplot_smooth 工具。
% 实现思路：FFT 后对幅度谱做分数倍频程几何平均，再以 dB 对数频率轴绘出。
%
% Author / 作者: Zhenxian Li

if nargin < 3 || isempty(fracOct)
    fracOct = 3;                      % default 1/3-oct / 默认 1/3 倍频程
end

x = x(:);
Nfft = 2^nextpow2(max(length(x), 1024));
X    = fft(x, Nfft);
mag  = abs(X(1:floor(Nfft/2)+1));
f    = (0:floor(Nfft/2))' * fs / Nfft;

% Convert to dB after a tiny floor to avoid log(0). / 加最小值避免 log(0)
magDb = 20*log10(mag + eps);

% Fractional-octave smoothing over a log-spaced ratio. / 对数比例下平滑
ratio  = 2^(1/(2*fracOct));
magSmoothed = magDb;
for k = 2:length(f)-1
    fLo = f(k) / ratio;
    fHi = f(k) * ratio;
    idx = (f >= fLo) & (f <= fHi);
    if any(idx)
        magSmoothed(k) = mean(magDb(idx));
    end
end

semilogx(f(2:end), magSmoothed(2:end), 'LineWidth', 1.4);
grid on; box on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([20 min(fs/2, 24000)]);
set(gca, 'XMinorGrid', 'on');
end
