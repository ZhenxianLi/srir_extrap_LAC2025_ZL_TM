function h = fig(varargin)
% FIG  Convenience wrapper for FIGURE used throughout the LAC scripts.
%      LAC 论文脚本里使用的 figure 快捷包装函数。
%
% Usage / 用法:
%   fig             % new white figure with default size / 新建白底图窗
%   fig(name)       % pass a name through to FIGURE / 透传 name 到 FIGURE
%   h = fig(...)    % return the handle / 返回图窗句柄
%
% The original LAC scripts use a custom `fig` shorthand from McKenzie's
% private collection. This drop-in replacement keeps the call sites working
% without dragging in that extra dependency.
% 原 LAC 脚本里大量使用 McKenzie 个人工具里的 fig 简写，这里给出等价替代，
% 避免引入额外依赖。
%
% Author / 作者: Zhenxian Li

h = figure(varargin{:});
set(h, 'Color', 'w');
% Use a reasonable default size; users can resize later.
% 给一个常用的默认尺寸，用户后续可自由调整。
pos = get(h, 'Position');
set(h, 'Position', [pos(1) pos(2) 720 360]);
end
