function status = installDependencies(varargin)
% INSTALLDEPENDENCIES  One-stop check (and optional auto-fetch) of every
%                      external resource this project needs.
%                      集中检查 / 按需下载本项目所有外部资源。
%
% The two main entry scripts (run_paper_figures, run_beamforming) call
% this once at start-up. They never deal with download URLs, install
% locations or third-party setup themselves — that lives only here.
% 两个主入口脚本（run_paper_figures、run_beamforming）只在最开头调用一次
% 本函数，下载链接、安装路径、第三方初始化全部封装在这里，不会散落各处。
%
% Resources covered / 涵盖的资源:
%   1) Git submodules under external/    — verified by sentinel files;
%                                          missing ones cause an error with
%                                          the exact git command to run.
%      external/ 下的 git 子模块            缺失时报错并打印需要运行的 git 命令。
%   2) Arni 6DoF SOFA dataset (~565 MB)  — auto-downloaded from the GitHub
%                                          Release mirror of the Zenodo
%                                          record (5720724).
%      Arni 6DoF SOFA 数据集（约 565 MB）   缺失时自动从 GitHub Release 镜像下载。
%   3) AMT toolbox (mckenzie2025)        — optional. Detected if already on
%                                          path; auto-downloaded + unzipped
%                                          only when 'InstallAMT' is true.
%      AMT 工具箱（mckenzie2025）           可选；已在 path 上就直接用，否则
%                                          仅当 InstallAMT=true 时下载安装。
%
% Name-Value options / 命名参数:
%   'InstallSofa'  (true)  — auto-download SOFA dataset when missing
%                            缺失时自动下载 SOFA 数据集
%   'InstallAMT'   (false) — auto-download + unzip AMT when missing
%                            缺失时自动下载并解压 AMT
%   'AMTRoot'      ('')    — target folder for AMT; default
%                            C:\MATLAB\amtoolbox (Windows) or
%                            ~/MATLAB/amtoolbox (Linux/macOS)
%                            AMT 安装目录（默认见上）
%   'Interactive'  (true)  — ask for confirmation before large downloads
%                            大文件下载前提示确认
%   'FailOnMissing' (true) — error out when a required resource is still
%                            missing after the auto-fetch step
%                            必需资源缺失时是否抛错（默认是）
%
% Output / 输出 (struct):
%   .submodules  logical — every required submodule is present
%                          所有必需子模块齐全
%   .sofaFile    char    — absolute path to the SOFA file (or '' if missing)
%                          SOFA 文件绝对路径（缺失则为 ''）
%   .amt         logical — mckenzie2025 is callable in this session
%                          mckenzie2025 已可用
%
% Examples / 示例:
%   installDependencies();                    % 一键检查 + 自动下 SOFA
%   installDependencies('InstallAMT', true);  % 顺便把 AMT 也装好
%   installDependencies('Interactive', false);% 静默（CI 用）
%
% Author / 作者: Zhenxian Li

% ----- parse options / 参数解析 -------------------------------------------
p = inputParser;
addParameter(p, 'InstallSofa',   true,  @islogical);
addParameter(p, 'InstallAMT',    false, @islogical);
addParameter(p, 'AMTRoot',       '',    @ischar);
addParameter(p, 'Interactive',   true,  @islogical);
addParameter(p, 'FailOnMissing', true,  @islogical);
parse(p, varargin{:});
opt = p.Results;

repoRoot = fileparts(fileparts(mfilename('fullpath')));     % utils/.. -> root
status   = struct('submodules', false, 'sofaFile', '', 'amt', false);

fprintf('\n========================================================\n');
fprintf(' Dependency check / 依赖检查\n');
fprintf('========================================================\n');

% ----- 1) submodules / 子模块 ---------------------------------------------
status.submodules = local_checkSubmodules(repoRoot);

% ----- 2) SOFA dataset / SOFA 数据集 --------------------------------------
status.sofaFile = local_ensureSofa(repoRoot, opt.InstallSofa);

% ----- 3) AMT (optional) / AMT（可选） -----------------------------------
status.amt = local_ensureAMT(opt.InstallAMT, opt.AMTRoot, opt.Interactive);

% ----- 4) Patch legacy SOFA conventions into SOFAtoolbox ------------------
% ----- 4) 把老版 SOFA convention patch 回 SOFAtoolbox --------------------
if status.submodules
    local_patchSofaConventions(repoRoot);
end

fprintf('========================================================\n\n');

% ----- enforce hard failures / 强制报错 -----------------------------------
if opt.FailOnMissing
    if ~status.submodules
        error('installDependencies:missingSubmodules', ...
            ['Required git submodules are missing. Run this in a shell:\n' ...
             '    git submodule update --init --recursive\n\n' ...
             '缺失必需子模块，请在 shell 中运行：\n' ...
             '    git submodule update --init --recursive']);
    end
    if opt.InstallSofa && isempty(status.sofaFile)
        error('installDependencies:sofaMissing', ...
            ['SOFA dataset still missing after auto-download attempt.\n' ...
             '自动下载后仍未找到 SOFA 数据集，请查看上方日志或手动下载。']);
    end
end
end


% =========================================================================
% local helpers / 子函数
% =========================================================================

function ok = local_checkSubmodules(repoRoot)
% Verify each submodule by looking for a known file inside it.
% 通过查找每个子模块里的标志文件来确认其已就位。
required = {
    'SOFAtoolbox',                       'SOFAstart.m';
    'Higher-Order-Ambisonics',           'getRSH.m';
    'Spherical-Harmonic-Transform',      'getSH.m';
    'Spherical-Array-Processing',        'sphPWDmap.m';
    'binaural_ambisonic_preprocessing',  'load_ambisonic_configuration.m'};

ok = true;
for i = 1:size(required, 1)
    name     = required{i, 1};
    sentinel = required{i, 2};

    % Prefer external/<name>; fall back to legacy bundled copy at root.
    % 优先 external/<name>，退路是仓库根目录下的同名文件夹（兼容老布局）。
    candidates = { fullfile(repoRoot, 'external', name), ...
                   fullfile(repoRoot, name) };

    found = false; where = '';
    for c = 1:numel(candidates)
        if isfolder(candidates{c}) && ...
                ~isempty(dir(fullfile(candidates{c}, '**', sentinel)))
            found = true; where = candidates{c}; break;
        end
    end

    if found
        fprintf('  [OK]   submodule  %-38s -> %s\n', name, where);
    else
        fprintf('  [MISS] submodule  %-38s  NOT FOUND\n', name);
        ok = false;
    end
end

if ~ok
    fprintf('\n  --> Run from a shell to pull submodules:\n');
    fprintf('      git submodule update --init --recursive\n');
    fprintf('  --> 请在 shell 中运行以拉取子模块：\n');
    fprintf('      git submodule update --init --recursive\n\n');
end
end


function local_patchSofaConventions(repoRoot)
% Copy legacy SOFA convention CSV files into the SOFAtoolbox submodule's
% conventions/ directory so SOFAstart can compile them to .mat. Needed
% because new SOFAtoolbox (2.x) dropped SingleRoomDRIR — but our 2021
% Arni dataset still declares that convention.
% 把老版 SOFA convention CSV 拷到 SOFAtoolbox submodule 的 conventions/，
% 让 SOFAstart 自动编译为 .mat。新版 SOFAtoolbox (2.x) 已经移除了
% SingleRoomDRIR，但 2021 年的 Arni 数据集仍然声明使用它。

srcDir = fullfile(repoRoot, 'utils', 'conventions');
if ~isfolder(srcDir); return; end                          % nothing to patch

dstDir = fullfile(repoRoot, 'external', 'SOFAtoolbox', 'SOFAtoolbox', 'conventions');
if ~isfolder(dstDir)
    dstDir = fullfile(repoRoot, 'SOFAtoolbox', 'conventions');     % legacy
end
if ~isfolder(dstDir)
    fprintf('  [SKIP] SOFAtoolbox conventions dir not found — convention patch skipped.\n');
    return;
end

csvFiles = dir(fullfile(srcDir, '*.csv'));
for k = 1:numel(csvFiles)
    name = csvFiles(k).name;
    src  = fullfile(srcDir, name);
    dst  = fullfile(dstDir, name);

    % Skip if already present and same size (cheap freshness check)
    % 已存在且大小相同就跳过（廉价的新鲜度检查）
    if exist(dst, 'file')
        s1 = dir(src); s2 = dir(dst);
        if s1.bytes == s2.bytes
            fprintf('  [OK]   SOFA patch  %-38s already in place\n', name);
            continue;
        end
    end

    try
        copyfile(src, dst, 'f');
        fprintf('  [PATCH] SOFA patch  %-38s -> conventions/\n', name);
        fprintf('          (SOFAstart will compile it to .mat next run / SOFAstart 下次启动自动编译为 .mat)\n');
    catch ME
        fprintf('  [WARN] SOFA patch  %-38s failed: %s\n', name, ME.message);
    end
end

% Tell the submodule's git to ignore the patched CSV locally — keeps
% `git status` clean without modifying the submodule's tracked .gitignore.
% 让 submodule 本地 git 忽略我们 patch 过去的 CSV，保持 git status 干净，
% 又不修改 submodule 自身被 tracked 的 .gitignore。
local_addToGitExclude(repoRoot, fullfile('external', 'SOFAtoolbox'), ...
    'SOFAtoolbox/conventions/SingleRoomDRIR_*.csv');
end


function local_addToGitExclude(repoRoot, submoduleRelPath, pattern)
% Resolve the submodule's real gitdir (it stores .git as a gitdir-file
% pointing into <super>/.git/modules/<path>) and append `pattern` to
% info/exclude if not already present. info/exclude is local-only and
% never committed.
% 解析 submodule 的真实 gitdir（submodule 里的 .git 是个 gitdir-file，
% 指向 <super>/.git/modules/<path>），然后把 pattern 加到 info/exclude
% （如果还没有的话）。info/exclude 是本地的，不会入库。

subRoot = fullfile(repoRoot, submoduleRelPath);
gitMeta = fullfile(subRoot, '.git');
gitDir  = '';
if exist(gitMeta, 'dir')
    gitDir = gitMeta;
elseif exist(gitMeta, 'file')
    raw = fileread(gitMeta);
    tok = regexp(strtrim(raw), '^gitdir:\s*(.+)$', 'tokens', 'once');
    if ~isempty(tok)
        candidate = strtrim(tok{1});
        if isfolder(fullfile(subRoot, candidate))
            gitDir = fullfile(subRoot, candidate);
        elseif isfolder(candidate)
            gitDir = candidate;
        end
    end
end
if isempty(gitDir); return; end

infoDir = fullfile(gitDir, 'info');
if ~isfolder(infoDir); mkdir(infoDir); end
excl = fullfile(infoDir, 'exclude');

existing = '';
if exist(excl, 'file'); existing = fileread(excl); end
if contains(existing, pattern); return; end

fid = fopen(excl, 'a');
if fid == -1; return; end
cleaner = onCleanup(@() fclose(fid));     %#ok<NASGU>
if ~isempty(existing) && existing(end) ~= newline
    fprintf(fid, '\n');
end
fprintf(fid, '# Added by installDependencies (patched SOFA convention)\n');
fprintf(fid, '%s\n', pattern);
end


function sofaFile = local_ensureSofa(repoRoot, autoFetch)
% Check the dataset file; if absent and autoFetch is true, hand off to
% downloadSrirDataset (the single source of truth for the URL).
% 检查数据集是否存在；如缺失且 autoFetch=true，则交给 downloadSrirDataset
% 下载（URL 集中放在那个函数里，避免散落）。
sofaName = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';
sofaPath = fullfile(repoRoot, '6dof_SRIRs_eigenmike_SH', sofaName);

if exist(sofaPath, 'file')
    sofaFile = sofaPath;
    fprintf('  [OK]   SOFA dataset present / 已存在\n         %s\n', sofaPath);
    return;
end

if ~autoFetch
    sofaFile = '';
    fprintf('  [MISS] SOFA dataset missing / 缺失\n');
    fprintf('         pass ''InstallSofa'', true to fetch automatically\n');
    fprintf('         如需自动下载请传入 ''InstallSofa'', true\n');
    return;
end

fprintf('  [DL]   SOFA dataset missing — downloading (~565 MB)\n');
fprintf('         SOFA 数据集缺失，开始下载（约 565 MB）...\n');
try
    sofaFile = downloadSrirDataset(repoRoot);
catch ME
    sofaFile = '';
    fprintf('  [FAIL] SOFA download failed / 下载失败: %s\n', ME.message);
end
end


function ok = local_ensureAMT(doInstall, amtRoot, interactive)
% Return true iff mckenzie2025 is callable from MATLAB right now.
% 当且仅当 mckenzie2025 当前 MATLAB 会话可调用时返回 true。

% --- already reachable? ---
if exist('mckenzie2025', 'file') == 2
    fprintf('  [OK]   AMT (mckenzie2025) on path / 已在 path 上\n');
    ok = true; return;
end

% --- pick install location ---
if isempty(amtRoot)
    if ispc
        amtRoot = 'C:\MATLAB\amtoolbox';
    else
        amtRoot = fullfile(getenv('HOME'), 'MATLAB', 'amtoolbox');
    end
end

% --- installed at amtRoot but not yet initialised? ---
% The zip may unpack into a versioned subfolder (e.g. amtoolbox-1.6.0/),
% so search for the real amt_start.m anywhere under amtRoot.
% zip 解压后通常会带一层版本号子目录（如 amtoolbox-1.6.0/），所以这里
% 在 amtRoot 下递归搜 amt_start.m，找到真正的 AMT 根目录。
amtStartDir = local_findAmtStart(amtRoot);
if ~isempty(amtStartDir)
    fprintf('  [INFO] AMT found at %s — running amt_start ...\n', amtStartDir);
    addpath(amtStartDir);
    try
        run(fullfile(amtStartDir, 'amt_start.m'));
    catch ME
        fprintf('  [WARN] amt_start raised an error: %s\n', ME.message);
    end
    ok = exist('mckenzie2025', 'file') == 2;
    if ok
        fprintf('  [OK]   AMT initialised / 已初始化\n');
    else
        fprintf('  [WARN] amt_start ran but mckenzie2025 still not visible.\n');
    end
    return;
end

% --- not installed; either skip or auto-install ---
if ~doInstall
    fprintf('  [SKIP] AMT (mckenzie2025) not installed — PBC-2 colouration plot will be skipped.\n');
    fprintf('         未安装 AMT —— 论文里 PBC-2 音色那张子图将被跳过。\n');
    fprintf('         To install automatically:  installDependencies(''InstallAMT'', true)\n');
    fprintf('         如需自动安装：             installDependencies(''InstallAMT'', true)\n');
    ok = false; return;
end

% --- auto-install path ---
% Folder on SourceForge is "AMT 1.x" (covers 1.5, 1.6, ...); zip carries
% the version in its filename. Auxdata is fetched on demand by amt_start,
% so the zip itself is small (~14 MB).
% SourceForge 上的文件夹叫 "AMT 1.x"（涵盖 1.5、1.6 等版本）；版本号体现
% 在 zip 文件名里。auxdata 由 amt_start 按需下载，所以 zip 本体只有 ~14 MB。
amtVersion = '1.6.0';
url = sprintf(['https://sourceforge.net/projects/amtoolbox/files/' ...
               'AMT%%201.x/amtoolbox-full-%s.zip/download'], amtVersion);
zipFile = fullfile(tempdir, sprintf('amtoolbox-full-%s.zip', amtVersion));

if interactive
    fprintf('\n  AMT %s will be downloaded (~14 MB; auxdata fetched on demand)\n', amtVersion);
    fprintf('  and unzipped to:\n    %s\n', amtRoot);
    fprintf('  即将下载 AMT %s（zip 约 14 MB；辅助数据按需下载）并解压到上面的目录。\n', amtVersion);
    resp = input('  Continue? / 继续？ [Y/n]: ', 's');
    if ~isempty(resp) && lower(resp(1)) ~= 'y'
        fprintf('  Aborted by user / 已取消。\n');
        ok = false; return;
    end
end

fprintf('  [DL]   Downloading AMT %s ...\n', amtVersion);
try
    websave(zipFile, url, weboptions('Timeout', 1800));
catch ME
    fprintf('  [FAIL] AMT download failed / 下载失败: %s\n', ME.message);
    fprintf('         Manual install: https://amtoolbox.org/download.php\n');
    ok = false; return;
end

if ~isfolder(amtRoot); mkdir(amtRoot); end
fprintf('  [UNZIP] Unpacking to %s ...\n', amtRoot);
try
    unzip(zipFile, amtRoot);
catch ME
    fprintf('  [FAIL] Unzip failed / 解压失败: %s\n', ME.message);
    ok = false; return;
end

% Locate the real AMT root (may be a versioned subfolder created by unzip).
% 解压后 amt_start.m 可能在版本号子目录里，递归搜一下。
amtStartDir = local_findAmtStart(amtRoot);
if isempty(amtStartDir)
    fprintf('  [FAIL] Could not find amt_start.m anywhere under %s\n', amtRoot);
    fprintf('         在 %s 下未找到 amt_start.m\n', amtRoot);
    ok = false; return;
end

addpath(amtStartDir);
try
    run(fullfile(amtStartDir, 'amt_start.m'));
catch ME
    fprintf('  [WARN] amt_start raised an error: %s\n', ME.message);
end

ok = exist('mckenzie2025', 'file') == 2;
if ok
    fprintf('  [OK]   AMT installed and initialised / 安装并初始化完成\n');
    fprintf('         AMT root: %s\n', amtStartDir);
else
    fprintf('  [WARN] AMT unpacked but mckenzie2025 still not callable — see warnings above.\n');
end
end


function amtStartDir = local_findAmtStart(searchRoot)
% Return the directory that directly contains amt_start.m, or '' if none.
% Tries the search root first, then a single recursive sweep.
% 优先看 searchRoot 是否直接含 amt_start.m；否则在其下递归搜。
amtStartDir = '';
if exist(fullfile(searchRoot, 'amt_start.m'), 'file')
    amtStartDir = searchRoot; return;
end
hits = dir(fullfile(searchRoot, '**', 'amt_start.m'));
if ~isempty(hits)
    amtStartDir = hits(1).folder;
end
end
