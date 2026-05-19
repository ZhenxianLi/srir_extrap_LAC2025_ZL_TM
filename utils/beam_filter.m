function srir_out = beam_filter(srir_in, azDeg, elDeg, beta, sigma)
% BEAM_FILTER  In-memory direction-selective filter for an HOA SRIR chunk.
%              对 HOA SRIR 短窗做"方向选择性过滤"的内存版本。
%
% This is the in-memory port of the core algorithm inside steerAmbiSOFA —
% identical SH-steering + Gaussian per-order taper + linear cross-fade —
% but operating directly on an Ns×C matrix so it can be called inside the
% LAC arrival loop without writing/reading SOFA files for every reflection.
% 本函数是 steerAmbiSOFA 核心算法的内存版本 —— 完全相同的 SH 波束 + 每阶
% 高斯锥窗 + 线性 cross-fade —— 但直接处理 Ns×C 矩阵，便于在 LAC 每个
% 反射的循环里逐窗调用（不必每次都写/读 SOFA 文件）。
%
% Usage / 用法:
%   srir_out = beam_filter(srir_in, azDeg, elDeg, beta, sigma)
%
% Inputs / 输入:
%   srir_in : Ns×C HOA matrix (ACN/N3D), one short arrival window
%             Ns×C 的 HOA 短窗（ACN/N3D 通道顺序）
%   azDeg   : look-direction azimuth in degrees   (0 = front, +90 = left, CCW)
%             波束指向方位角（度，0=前，+90=左）
%   elDeg   : look-direction elevation in degrees (+90 = zenith)
%             波束指向仰角（度，+90=顶）
%   beta    : cross-fade ratio, 0…1
%             0 = pass-through original chunk
%             1 = pure beam-only output
%             cross-fade 比例：0=完全原始，1=完全波束
%   sigma   : Gaussian beam-width parameter (>0)
%             σ ≈ 1   → narrow beam
%             σ ≈ 2–3 → wider, softer beam
%             高斯波束宽度（σ 越小波束越窄）
%
% Output / 输出:
%   srir_out : Ns×C filtered HOA matrix (same size as input)
%              过滤后的 Ns×C HOA 矩阵（与输入同尺寸）
%
% Algorithm / 算法（与 steerAmbiSOFA 完全一致）:
%   1) Build a real ACN/N3D SH steering vector Y at (az, el).
%      在 (az, el) 处构造实数 ACN/N3D 球谐 steering 向量 Y。
%   2) Apply per-order Gaussian taper w_n = exp(-n/σ) to widen/narrow beam.
%      用每阶高斯权重 w_n = exp(-n/σ) 控制波束宽度。
%   3) Decode to plane-wave: s(t) = X · Ywᵀ        (Ns×1)
%      解码到一路平面波信号
%   4) Re-encode to HOA:     X_beam = s(t) · Yw    (Ns×C)
%      重新编码回 HOA
%   5) Cross-fade: out = (1-β)·X + β·X_beam
%      线性混合：(1-β)·原始 + β·波束
%
% Author / 作者: Zhenxian Li

[~, C] = size(srir_in);
N      = round(sqrt(C) - 1);                       % HOA order, since C = (N+1)^2

% --- 1) Real ACN/N3D SH steering vector at (az, el) ----------------------
Y  = local_getRealSH_N3D(N, deg2rad(azDeg), deg2rad(elDeg));    % 1×C

% --- 2) Per-order Gaussian taper, expanded to channel-vector --------------
w_per_order = exp( -(0:N).' / sigma );                          % (N+1)×1
w_vec       = local_expandWeights_ACN(w_per_order);             % 1×C
Yw          = Y .* w_vec;
Yw          = Yw / norm(Yw);                                    % unit energy

% --- 3) Decode -> single plane wave ---------------------------------------
% --- 4) Re-encode -> HOA --------------------------------------------------
s          = srir_in * Yw.';                                    % Ns×1
srir_beam  = s * Yw;                                            % Ns×C

% --- 5) Linear cross-fade -------------------------------------------------
srir_out   = (1 - beta) * srir_in + beta * srir_beam;
end


% =========================================================================
% local helpers — vendored from steerAmbiSOFA.m to keep beam_filter
% self-contained. 与 steerAmbiSOFA.m 完全相同的两个辅助函数。
% =========================================================================

function Y = local_getRealSH_N3D(N, az, el)
% Real, N3D-normalized SH coefficients in ACN order at (az, el).
% Returns a 1×(N+1)^2 unit-energy row vector.
theta = pi/2 - el;
nCh   = (N+1)^2;
Y     = zeros(1, nCh);

for n = 0:N
    P = legendre(n, cos(theta), 'sch');    % Schmidt semi-normalised
    for m = -n:n
        idx = n*(n+1) + m + 1;             % ACN channel index
        Pnm = P(abs(m)+1);
        K   = sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
        if m > 0
            Y(idx) = sqrt(2) * K * Pnm * cos( m * az );
        elseif m < 0
            Y(idx) = sqrt(2) * K * Pnm * sin( abs(m) * az );
        else
            Y(idx) =         K * Pnm;
        end
    end
end
Y = Y / norm(Y);
end


function wvec = local_expandWeights_ACN(w_order)
% Convert (N+1)×1 per-order weights to a 1×(N+1)^2 channel-aligned vector.
N    = numel(w_order) - 1;
wvec = zeros(1, (N+1)^2);
idx  = 1;
for n = 0:N
    numM = 2*n + 1;
    wvec(idx : idx+numM-1) = w_order(n+1);
    idx = idx + numM;
end
end
