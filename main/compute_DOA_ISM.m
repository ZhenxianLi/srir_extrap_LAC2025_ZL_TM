
%subspace Decomposition
%cut direct sound
%scale direct sound
%rotate direct sound
%rejoin




% from single recording
%--------------------------------------------------------------------------------------------------------------------------
%% load and plot sofa file

% 包括对早期反射进行缩放并通过改变直达声到达时间与幅度来模拟在房间内
% 不同位置处的空间音频响应。

clc;
clear all;
close all;

% Start SOFA and load function library
% 加载SOFA相关API以及其他所需函数库，用于处理高阶Ambisonics和SOFA格式的脉冲响应

% get the parentFolder path
% 获取当前工作目录所在的父文件夹路径，便于添加需要的依赖路径
parentFolder = fileparts(pwd);

addpath(fullfile(parentFolder, 'API_MO'));
addpath(fullfile(parentFolder, 'Spherical-Harmonic-Transform-master'));
addpath(fullfile(parentFolder, 'Higher-Order-Ambisonics-master'));
addpath(genpath(fullfile(parentFolder, 'binaural_ambisonic_preprocessing-main')));
addpath(fullfile(parentFolder, 'audio_samples'));
addpath(fullfile(parentFolder, 'SRIR-Subspace-Decomposition-master'));
SOFAstart;


%run this to load decoder
% 如果之前没有加载Ambisonics配置，则运行load_ambisonic_configuration来获得Ambisonics解码器
if ~exist('aio_flag')
    load_ambisonic_configuration
end

% Load Variable Acoustics 6DoF Dataset (most reverberant set)
% 这里加载六自由度的SRIRs数据（最混响的版本），文件保存在 6dof_SRIRs_eigenmike_SH 文件夹中

parentFolder = fileparts(pwd);
irPath1 = fullfile(parentFolder, '6dof_SRIRs_eigenmike_SH/');
irName1 = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';

% download the SRIR.sofa from Github Release
% 如果本地没有该SOFA文件则从Github下载
url = 'https://github.com/ZhenxianLi/ZHENXIAN_SRIRextrapolation/releases/download/6dof_SRIRs_eigenmike_SH/6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';
% Construct the full path for the saved file
outputFile = fullfile(irPath1, irName1);
% Check if the file already exists to avoid duplicate downloads
if ~exist(outputFile, 'file')
    fprintf('start Download the file: %s\n', irName1);
    try
        % download to pass
        websave(outputFile, url);
        fprintf('downloadSucess: %s\n', outputFile);
    catch ME
        fprintf('downloadFail: %s\n', ME.message);
        return; % 终止执行
    end
else
    fprintf('File already exists, skipping download: %s\n', outputFile);
end

% 使用SOFAAPI加载该SOFA文件
sofa1 = SOFAload([irPath1,irName1]);
fs = sofa1.Data.SamplingRate;
%SOFAplotGeometry(sofa1);   % plot the source and listen position
close all;
%% check inf
% 这里查看SOFA文件中的一些元数据和尺寸信息
sofa1.GLOBAL_RoomType;
sofa1.GLOBAL_SOFAConventionsVersion;
sofa1.ListenerPosition;
sofa1.SourcePosition;
sofa1.EmitterPosition;
size(sofa1.Data);


%%

% 示例：定义房间、声源、麦克风位置，声速 = 343 m/s
roomSize = [7.87, 5.75, 2.91];       % 房间尺寸 (长, 宽, 高)
sourcePos = sofa1.SourcePosition(13,:);  % 声源坐标
micPos    = sofa1.ListenerPosition(13,:);      % 麦克风坐标
maxRef    = 1;             % 只计算 maxRef 次镜像
c         = 343;           % 声速

% 调用函数
T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxRef, c);
T.Arr
% 查看输出表格 T
% 包含 Nx,Ny,Nz, Ximg,Yimg,Zimg, Distance, ArrTime, Theta, Phi



function T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxReflections, c)
% COMPUTE_DOA_ISM_WITH_PLOT  使用镜像源法计算并可视化 DOA、到达时间等信息
%
% T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, maxReflections, c)
%
% 输入参数：
%   roomSize       : [Lx, Ly, Lz] 房间尺寸 (米)
%   sourcePos      : [xs, ys, zs] 声源坐标 (米)
%   micPos         : [xm, ym, zm] 麦克风坐标 (米)
%   maxReflections : 最大镜像次数 (推荐 1~2)
%   c              : 声速 (默认 343 m/s)
%
% 输出：
%   T : 表格，包含以下列：
%       Nx, Ny, Nz : 镜像次数(分别代表 x,y,z 方向的镜像index)
%       Ximg, Yimg, Zimg : 镜像源位置
%       Distance : 从镜像源到麦克风的距离 (米)
%       ArrTime  : 到达时间 (秒)
%       Theta    : 方位角 (Azimuth, 弧度)
%       Phi      : 仰角 (Elevation, 弧度)
%
% 演示：
%   roomSize = [5, 4, 3];
%   sourcePos = [2, 1.5, 1.2];
%   micPos    = [3, 2, 1];
%   T = compute_DOA_ISM_with_plot(roomSize, sourcePos, micPos, 1, 343);

    if nargin < 5
        c = 343;  % 如果未指定声速，默认 343 m/s
    end
    
    % 解析房间尺寸
    Lx = roomSize(1);
    Ly = roomSize(2);
    Lz = roomSize(3);

    xs = sourcePos(1); ys = sourcePos(2); zs = sourcePos(3);
    xm = micPos(1);   ym = micPos(2);   zm = micPos(3);

    % 预分配结果
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

    % 计数器：一次反射、二次反射
    singleCount = 0;  
    doubleCount = 0;  
    
    figure('Name','ISM Visualization','NumberTitle','off'); 
    hold on;
    grid on;
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Image Source Method (ISM) Visualization');
    
    % 绘制房间边界 (可选) - 简单画一个线框
    drawRoomFrame(0,0,0, Lx,Ly,Lz);

    % 绘制麦克风位置
    plot3(xm, ym, zm, 'ks', 'MarkerFaceColor','k','MarkerSize',8);
    text(xm, ym, zm, '  Mic', 'Color','k','FontWeight','bold');
    
    % 绘制原始声源 (nx=0,ny=0,nz=0)
    plot3(xs, ys, zs, 'r^', 'MarkerFaceColor','r','MarkerSize',8);
    text(xs, ys, zs, '  Source', 'Color','r','FontWeight','bold');

    % 遍历镜像次数
    for nx = -maxReflections:maxReflections
        for ny = -maxReflections:maxReflections
            for nz = -maxReflections:maxReflections

                % sumOrder = abs(nx)+abs(ny)+abs(nz)
                sumOrder = abs(nx) + abs(ny) + abs(nz);

                % 如果 sumOrder 在 [1, maxReflections]，则表示一次或二次(或更多)反射
                % 这里 maxReflections=2 时, sumOrder=1 ->一次, sumOrder=2->二次
                if sumOrder >= 1 && sumOrder <= maxReflections
                    % 计算镜像源位置
                    x_img = xs + nx * Lx;
                    y_img = ys + ny * Ly;
                    z_img = zs + nz * Lz;
                
                    % 方向向量：镜像源 -> 麦克风
                    d = [x_img - xm, y_img - ym, z_img - zm];
                
                    % 计算距离
                    dist = norm(d); % 欧几里得距离
                    % 计算到达时间
                    arrTime = dist / c;
                
                    % 计算 DOA (方位角 & 仰角)
                    theta = atan2(d(2), d(1));  % 方位角
                    phi   = atan2(d(3), sqrt(d(1)^2 + d(2)^2)); % 仰角

                    % 存储结果
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
                
                    % 区分一次反射( sumOrder=1 )与二次反射( sumOrder=2 )
                    if sumOrder == 1
                        singleCount = singleCount + 1;
                        % 绘制成绿色圆点(示例)
                        plot3(x_img, y_img, z_img, 'go','MarkerFaceColor','g','MarkerSize',5);
                    elseif sumOrder == 2
                        doubleCount = doubleCount + 1;
                        % 绘制成蓝色圆点(示例)
                        plot3(x_img, y_img, z_img, 'bo','MarkerFaceColor','b','MarkerSize',5);
                    end

                    % 连线：麦克风 -> 镜像源
                    line([xm x_img],[ym y_img],[zm z_img], ...
                         'Color',[0.3 0.3 0.8],'LineStyle','--');
                end
            end
        end
    end
    
    % 输出表格
    T = table(Nx_arr, Ny_arr, Nz_arr, ...
              Ximg_arr, Yimg_arr, Zimg_arr, ...
              Dist_arr, ArrTime_arr, Theta_arr, Phi_arr, ...
        'VariableNames', {'Nx','Ny','Nz','Ximg','Yimg','Zimg','Distance','ArrTime','Theta','Phi'});
    
    % 在 Command Window 中显示表格
    disp(T);

    % 显示统计结果
    fprintf('\n--- Reflection Stats ---\n');
    fprintf('  一次反射数量: %d\n', singleCount);
    fprintf('  二次反射数量: %d\n', doubleCount);
    
    % 调整视角
    view(3);
    hold off;
end

% 辅助函数：绘制房间外框
function drawRoomFrame(x0, y0, z0, Lx, Ly, Lz)
% 在当前坐标轴绘制一个线框，用于简单展示房间
    X = [x0, x0+Lx, x0+Lx, x0,     x0,     x0+Lx, x0+Lx, x0];
    Y = [y0, y0,     y0+Ly, y0+Ly, y0,     y0,     y0+Ly, y0+Ly];
    Z = [z0, z0,     z0,     z0,     z0+Lz, z0+Lz, z0+Lz, z0+Lz];
    edges = [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4];
    for i=1:length(edges)-1
       plot3([X(edges(i)) X(edges(i+1))], [Y(edges(i)) Y(edges(i+1))], ...
             [Z(edges(i)) Z(edges(i+1))], 'k-','LineWidth',1.0);
       hold on;
    end
end