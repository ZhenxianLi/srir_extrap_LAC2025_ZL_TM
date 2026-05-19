function repoRoot = setupPaths()
% SETUPPATHS  Add every dependency on the MATLAB path for this project.
%             加入本项目所有依赖到 MATLAB path。
%
% Output / 输出:
%   repoRoot - absolute path of the repository root / 仓库根目录绝对路径
%
% Lookup order for each dependency / 每个依赖的查找顺序:
%   1) external/<name>           (preferred — populated via git submodule)
%                                 推荐：通过 git submodule 拉取
%   2) <name> at repository root (legacy bundled copy, for backwards compat)
%                                 兼容老布局：依赖直接在仓库根目录里
%
% After this you can call `SOFAstart`, `rotateHOA_N3D`, `getRSH`, etc.
% Also adds the project's own utils/, FinalCode/, and audio_samples/ folders.
% 调用后即可直接使用 SOFAstart / rotateHOA_N3D / getRSH 等函数，
% 同时也把项目自带的 utils/、FinalCode/、audio_samples/ 加入 path。
%
% Author / 作者: Zhenxian Li

repoRoot = fileparts(fileparts(mfilename('fullpath')));   % utils/.. -> repo root

% --- project's own folders / 项目自带文件夹 ---------------------------------
addpath(fullfile(repoRoot, 'utils'));
addpath(fullfile(repoRoot, 'FinalCode'));
addpath(fullfile(repoRoot, 'audio_samples'));

% --- external dependencies / 外部依赖 ---------------------------------------
deps = struct( ...
    'sofaToolbox',  {{'SOFAtoolbox', 'API_MO'}}, ...
    'hoaLib',       {{'Higher-Order-Ambisonics', 'Higher-Order-Ambisonics-master'}}, ...
    'shtLib',       {{'Spherical-Harmonic-Transform', 'Spherical-Harmonic-Transform-master'}}, ...
    'sapLib',       {{'Spherical-Array-Processing', 'Spherical-Array-Processing-master'}}, ...
    'binPrep',      {{'binaural_ambisonic_preprocessing', 'binaural_ambisonic_preprocessing-main'}});

fields = fieldnames(deps);
for i = 1:numel(fields)
    candidates = deps.(fields{i});
    found = false;
    for c = 1:numel(candidates)
        % prefer external/<name> (submodule); fall back to <name> at root
        % 优先 external/<name>（submodule），退路是仓库根目录的 <name>
        p_ext   = fullfile(repoRoot, 'external', candidates{c});
        p_root  = fullfile(repoRoot, candidates{c});
        if isfolder(p_ext)
            addPathRecursive(p_ext, fields{i});
            found = true; break;
        elseif isfolder(p_root)
            addPathRecursive(p_root, fields{i});
            found = true; break;
        end
    end
    if ~found
        warning('setupPaths:missingDep', ...
            ['Dependency "%s" not found. Did you run ' ...
             '`git submodule update --init --recursive`?\n' ...
             '未找到依赖 "%s"，请先执行 git submodule update --init --recursive'], ...
             candidates{1}, candidates{1});
    end
end

% Special case: SOFAtoolbox layout differs slightly — the runtime sits in
% SOFAtoolbox/. The newer toolbox keeps the same entry name `SOFAstart`.
% 特殊处理：SOFAtoolbox 的目录结构略不同，但入口仍是 SOFAstart。

end

function addPathRecursive(p, depName)
% Recursive addpath but skip vendored test scripts that may clear globals.
% 递归 addpath，但跳过会 clear globals 的 vendored test 脚本。
allDirs = genpath(p);
addpath(allDirs);
fprintf('  added (%s): %s\n', depName, p);
end
