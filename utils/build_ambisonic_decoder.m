function build_ambisonic_decoder(repoRoot)
% BUILD_AMBISONIC_DECODER  Build SH_ambisonic_binaural_decoder in the base
%                          workspace, robust to the bundled diagnostics
%                          tail that lives in the submodule.
%                          在 base workspace 里构建双耳 Ambisonic 解码器，
%                          并对 submodule 自带的诊断尾巴做容错处理。
%
% Why a wrapper / 为什么要包一层：
%   `binaural_ambisonic_preprocessing/test_scripts/load_ambisonic_configuration.m`
%   ends (line 103) with an unconditional call to `test_ambisonic_decoder` —
%   a *diagnostics* script that:
%     1. uses a relative path `dir('hrirs/*.wav')` (needs cwd == submodule root);
%     2. depends on `mckenzie2022` from AMT, which is optional here.
%   The decoder itself is fully built before that tail runs, so we
%   cd into the submodule root and swallow any failure of the tail.
%
%   该 submodule 自带的 `load_ambisonic_configuration.m` 在第 103 行
%   无条件调用诊断脚本 `test_ambisonic_decoder`，后者：
%     1) 用相对路径 dir('hrirs/*.wav') 加载数据（依赖 cwd 是 submodule 根目录）；
%     2) 依赖 AMT 的 mckenzie2022 模型（本项目里是可选依赖）。
%   而真正的 decoder 在诊断尾巴跑之前就已经构建好，因此本函数：
%     - cd 到 submodule 根目录，让相对路径有效；
%     - 用 try/catch 吞掉诊断尾巴可能的错误，主流程不受影响。
%
% Input / 输入:
%   repoRoot - absolute path to the repository root / 仓库根目录绝对路径
%
% Output / 输出:
%   none — leaves SH_ambisonic_binaural_decoder (and the other variables
%   load_ambisonic_configuration creates) in the base workspace, ready for
%   the LAC reproduction script to consume.
%   无返回值；解码器及相关变量留在 base workspace，供论文复现脚本使用。
%
% Author / 作者: Zhenxian Li

prevWD = pwd;
cdGuard = onCleanup(@() cd(prevWD));     %#ok<NASGU>  % restore cwd on exit

submoduleRoot = fullfile(repoRoot, 'external', 'binaural_ambisonic_preprocessing');
if ~isfolder(submoduleRoot)
    error('build_ambisonic_decoder:missingSubmodule', ...
        ['binaural_ambisonic_preprocessing submodule not found at:\n  %s\n' ...
         'Run `git submodule update --init --recursive` first.\n' ...
         '未找到 binaural_ambisonic_preprocessing 子模块，请先运行 git submodule update。'], ...
        submoduleRoot);
end
cd(submoduleRoot);

% Snapshot base workspace because load_ambisonic_configuration line 41
% runs `clear variables`, which would wipe the caller's repoRoot, figDir,
% etc. We save now and reload after the script returns, so the caller's
% state is intact and newly created variables (decoder) survive too.
% load_ambisonic_configuration 第 41 行有 `clear variables`，会把调用者
% （脚本）的 repoRoot / figDir 等都清掉。先把 base 快照存盘，跑完再
% 重新 load 回去 —— 调用者的变量回来了，新生成的 decoder 也保留。
snapshotFile = [tempname() '.mat'];
restoreGuard = onCleanup(@() iSafeDelete(snapshotFile));    %#ok<NASGU>
hasAny = evalin('base', '~isempty(who)');
if hasAny
    assignin('base', 'i_buildDecoder_snap_', snapshotFile);
    evalin('base', 'save(i_buildDecoder_snap_); clear i_buildDecoder_snap_;');
end

try
    evalin('base', 'load_ambisonic_configuration; aio_flag = 1;');
catch ME
    % Decoder may have been built before the diagnostics tail crashed.
    % If it's already in base, the failure is benign — warn and continue.
    % 诊断尾巴失败前，decoder 通常已经构建好；只要它在 base 里，就视为良性失败。
    decoderInBase = evalin('base', ...
        'exist(''SH_ambisonic_binaural_decoder'', ''var'') == 1');
    if decoderInBase
        warning('buildDecoder:diagnosticsSkipped', ...
            ['Decoder built OK. Bundled diagnostics script\n' ...
             '(test_ambisonic_decoder) failed — safe to ignore.\n' ...
             'Decoder 已成功构建，自带诊断脚本失败可以忽略。\n' ...
             'Reason / 原因: %s'], ME.message);
    else
        if hasAny
            assignin('base', 'i_buildDecoder_snap_', snapshotFile);
            evalin('base', 'load(i_buildDecoder_snap_); clear i_buildDecoder_snap_;');
        end
        rethrow(ME);
    end
end

if hasAny
    % load preserves variables that already exist in base if they were
    % present in the .mat file; the decoder (new variable) stays put.
    % load 会把 .mat 里保存过的旧变量都还原回 base；新生成的 decoder
    % 不在 .mat 里，因此 load 不会动它。
    assignin('base', 'i_buildDecoder_snap_', snapshotFile);
    evalin('base', 'load(i_buildDecoder_snap_); clear i_buildDecoder_snap_;');
end

assert(evalin('base', 'exist(''SH_ambisonic_binaural_decoder'', ''var'') == 1'), ...
    'build_ambisonic_decoder:notBuilt', ...
    'SH_ambisonic_binaural_decoder was not created in base workspace.');
end


function iSafeDelete(p)
if exist(p, 'file'); delete(p); end
end
