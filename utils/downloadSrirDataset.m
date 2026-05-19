function sofaFile = downloadSrirDataset(repoRoot)
% DOWNLOADSRIRDATASET  Make sure the 6DoF Arni 50%-absorber SOFA file is present.
%                      确保 6DoF Arni 50% 吸声 SOFA 数据集已就位。
%
% Usage / 用法:
%   sofaFile = downloadSrirDataset();          % uses current repo
%   sofaFile = downloadSrirDataset(repoRoot);  % explicit root
%
% On first run this downloads the ~565 MB dataset from the project's
% GitHub Release mirror (which fronts the original Zenodo record:
%   https://zenodo.org/records/5720724 ).
% 第一次运行会从作者的 GitHub Release 镜像下载约 565 MB 的数据集
% （上游来源为 Zenodo: https://zenodo.org/records/5720724 ）。
%
% Author / 作者: Zhenxian Li

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

irDir   = fullfile(repoRoot, '6dof_SRIRs_eigenmike_SH');
irName  = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';
sofaFile = fullfile(irDir, irName);

if exist(sofaFile, 'file')
    fprintf('[downloadSrirDataset] OK already present / 已存在: %s\n', sofaFile);
    return;
end

if ~isfolder(irDir); mkdir(irDir); end

url = ['https://github.com/ZhenxianLi/ZHENXIAN_SRIRextrapolation/' ...
       'releases/download/6dof_SRIRs_eigenmike_SH/' irName];

fprintf('[downloadSrirDataset] downloading ~565 MB / 正在下载（约 565 MB）...\n');
fprintf('  src: %s\n  dst: %s\n', url, sofaFile);

try
    opts = weboptions('Timeout', 600);
    websave(sofaFile, url, opts);
    fprintf('[downloadSrirDataset] OK done / 完成: %s\n', sofaFile);
catch ME
    error('downloadSrirDataset:failed', ...
        ['Download failed (%s).\n' ...
         '下载失败 (%s)。\n' ...
         'You can place the file manually at:\n' ...
         '可手动下载并放到：\n  %s\n' ...
         'Original source / 原始来源: https://zenodo.org/records/5720724'], ...
        ME.message, ME.message, sofaFile);
end
end
