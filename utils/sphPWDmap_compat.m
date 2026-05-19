function [P_pwd, est_dirs, est_dirs_P] = sphPWDmap_compat(sphCOV, grid_dirs, nSrc, kappa)
% SPHPWDMAP_COMPAT  Backwards-compatibility shim for the legacy 4-input,
%                   3-output sphPWDmap signature.
%                   老版本 sphPWDmap（4 入 3 出）签名的兼容包装。
%
% Upstream sphPWDmap (Spherical-Array-Processing submodule, polarch) was
% simplified to 3 inputs / 2 outputs in newer commits — `kappa` is now
% hard-coded to 20 internally, and the per-peak power vector `est_dirs_P`
% is no longer returned. This project's paper code was written against
% the older API, so we re-implement the dropped pieces here:
%   * call the new upstream sphPWDmap to get the (unchanged) SRP map;
%   * redo peak picking ourselves to honour the caller's `kappa` AND to
%     return the power at each peak (est_dirs_P).
% Signature matches the legacy version exactly.
%
% 上游 sphPWDmap 在新版被简化成「3 入 2 出」——内部把 kappa 写死成 20，
% 也不再返回每个峰值方向的功率 est_dirs_P。论文代码按老 API 写的，所以
% 这里包一层：用新版 sphPWDmap 算 SRP 功率图，然后自己重新做峰值挑选，
% 既尊重调用方传入的 kappa，又能返回 est_dirs_P。签名与老版完全一致。
%
% Inputs / 输入:
%   sphCOV    : (order+1)^2 x (order+1)^2 SH-domain covariance matrix
%   grid_dirs : K x 2 [azi elev] in radians (typically from grid2dirs)
%   nSrc      : number of peaks to pick / 要找的峰值个数
%   kappa     : Von-Mises concentration; controls how aggressively previous
%               peaks are suppressed before picking the next one.
%               Von-Mises 浓度参数；越大每个峰挖洞越窄。
%
% Outputs / 输出:
%   P_pwd      : K x 1 SRP power map over grid_dirs
%   est_dirs   : nSrc x 2 [azi elev] of estimated peak directions
%   est_dirs_P : nSrc x 1 power at each estimated peak
%
% Author / 作者: Zhenxian Li

if nargin < 4 || isempty(kappa); kappa = 20; end

% 1) Power map — unchanged from upstream.
% 1) 功率图 —— 直接调用上游。
[P_pwd, ~] = sphPWDmap(sphCOV, grid_dirs, nSrc);

% 2) Re-do peak picking with the caller's kappa and capture per-peak powers.
% 2) 用调用方的 kappa 重做峰值挑选，并记录每个峰的功率。
grid_xyz     = unitSph2cart(grid_dirs);
P_residual   = P_pwd;
est_dirs     = zeros(nSrc, 2);
est_dirs_P   = zeros(nSrc, 1);
for k = 1:nSrc
    [~, peak_idx] = max(P_residual);
    est_dirs(k, :) = grid_dirs(peak_idx, :);
    est_dirs_P(k)  = P_pwd(peak_idx);

    VM_mean = grid_xyz(peak_idx, :);
    VM_mask = kappa / (2*pi*exp(kappa) - exp(-kappa)) * ...
              exp(kappa * grid_xyz * VM_mean');
    VM_mask = 1 ./ (1e-5 + VM_mask);
    P_residual = P_residual .* VM_mask;
end
end
