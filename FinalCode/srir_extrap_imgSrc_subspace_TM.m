% Spatial room impulse response extrapolation using image source reflections and
% subspace decomposition from single spatial room impulse response
%
% Zhenxian Li and Thomas McKenzie, 2025
% INSA de Lyon, France and University of Edinburgh, UK


% clc;
clear all;
% close all;

%% Settings

save_figs = true;
dim_room = [7.87 5.75 2.91]; % from the ArXiv paper
ord_imgsrc = 2; % what order of image source to compute

% Which SRIRs as original and target:
% % original
% idx_LS_src_orig=1;
% idx_LS_rec_orig=1;
% % target
% idx_LS_src_tar=1;
% idx_LS_rec_tar=1;

% % original
idx_LS_src_orig=2;
idx_LS_rec_orig=4;
% target
idx_LS_src_tar=3;
idx_LS_rec_tar=2;

% num samples before calculated start sample of arrival (for arrival extraction)
delay_arrival_initial_samp = 10;
arrival_dur_ms = 1; % how long to extract in ms (could do longer window for DS and shorter for others)
delay_arrival_initial_samp_DS = 50;
arrival_dur_ms_DS = 4; % how long to extract in ms (could do longer window for DS and shorter for others)

% Flags
flag_render_binaural = false;
flag_subspace_decomp = true;

%% Init
% Start SOFA and load function library
% get the new path
parentFolder = fileparts(pwd);

addpath(fullfile(parentFolder, 'API_MO'));
addpath(fullfile(parentFolder, 'Spherical-Harmonic-Transform-master'));
addpath(fullfile(parentFolder, 'Higher-Order-Ambisonics-master'));
addpath(genpath(fullfile(parentFolder, 'binaural_ambisonic_preprocessing-main')));
addpath(fullfile(parentFolder, 'audio_samples'));
% SOFAstart;

load_ambisonic_configuration

% Load Variable Acoustics 6DoF Dataset (medium reverberant set)
parentFolder = fileparts(pwd);
irPath1 = fullfile(parentFolder, '6dof_SRIRs_eigenmike_SH/');
irName1 = '6DoF_SRIRs_eigenmike_SH_50percent_absorbers_enabled.sofa';

%download the SRIR.sofa from Github Release
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

sofa1 = SOFAload([irPath1,irName1]);
fs = sofa1.Data.SamplingRate;





%% check inf
% sofa1.GLOBAL_RoomType;
% sofa1.GLOBAL_SOFAConventionsVersion;
% sofa1.ListenerPosition;
% sofa1.SourcePosition;
% sofa1.EmitterPosition;
% size(sofa1.Data);

% %check by hand
% firstReflectionIndex=[487;697;674;528;579;621;606;515;616;726;499;658;861;543;742;939;971;1146;1071;900;1171];

%% Get original and target SRIRs

idx_LS_srcrec_orig=idx_LS_rec_orig*3-3+idx_LS_src_orig; % original
idx_LS_srcrec_tar=idx_LS_rec_tar*3-3+idx_LS_src_tar; % target

%read IR and postion
srir_orig= squeeze(sofa1.Data.IR(idx_LS_srcrec_orig,:,:))'; % original IR
srir_tar= squeeze(sofa1.Data.IR(idx_LS_srcrec_tar,:,:))'; % target IR (for testing!)

% get source and receiver positions
pos_src_orig=sofa1.SourcePosition(idx_LS_srcrec_orig,:);
pos_list_orig=sofa1.ListenerPosition(idx_LS_srcrec_orig,:);

pos_src_tar=sofa1.SourcePosition(idx_LS_srcrec_tar,:);
pos_list_tar=sofa1.ListenerPosition(idx_LS_srcrec_tar,:);

% Plot geometry and label the original and target positions
SOFAplotGeometry(sofa1);   % plot the source and listen position

% text(pos_src_orig(1),pos_src_orig(2),'Orig','Color','red')
% text(pos_list_orig(1),pos_list_orig(2),'Orig','Color','red')
% text(pos_src_tar(1),pos_src_tar(2),'Tar','Color','green')
% text(pos_list_tar(1),pos_list_tar(2),'Tar','Color','green')

line([pos_src_orig(1) pos_list_orig(1)],[pos_src_orig(2) pos_list_orig(2)],'Color','blue')
line([pos_src_tar(1) pos_list_tar(1)],[pos_src_tar(2) pos_list_tar(2)],'Color','red')
% legend({'Listener Pos','Listener View Direction','Source Pos','Source View Direction','Original SRIR','Target SRIR'})
plot_time_vector = 0:1/fs:(size(srir_orig,1)/fs-1/fs);


%% intensity vector to detect source direction for original (measured)
% this will help to see if the image sources need rotation in any way
% (e.g. 90 degree rotation)

% find time of direct sound (in samples) - of measured SRIR
[delay_DS_orig_meas,~,~,del_orig_meas] = findDirectSound(srir_orig(:,1));

ds_duration = 0.004; % 4ms
ds_start_samp = delay_DS_orig_meas - 0.001*fs;
ds_end_samp = ds_start_samp + ds_duration*fs;

% intensity vector DoA analysis --- get measured direct sound direction 
dir_DS_meas_orig = doa_iv(srir_orig(ds_start_samp:ds_end_samp,:)); % in degs

%% SUBSPACE DECOMPOSITION
% separate out the direct and early reflections from the residual / diffuse!

if flag_subspace_decomp
    % parameters for the subspace decomposition, see the function header of srirSubspaceDecomp for details
    blockLenSmp = 32;
    hopSizeSmp = blockLenSmp / 8;
    kappa = 2.5;
    numBlocksGsvSumAvg = 32;
    residualEstimateLengthMs = 20;
    decompositionTimeLimitMs = 100;
    numBlocksSmoothThresh = 1;

    % samples x channels
    [srir_direct, srir_resid, numDirSubspaceComponents, gsvs, detectionThreshold, gsvSum, avgGsvSum] = ...
        srirSubspaceDecomp(srir_orig, fs, blockLenSmp, hopSizeSmp, kappa, numBlocksGsvSumAvg, residualEstimateLengthMs, ...
        decompositionTimeLimitMs, numBlocksSmoothThresh);
end

%% Image source detection

% original
[dist_imgsrc_orig, delay_imgsrc_orig, dir_imgsrc_orig, dist_DS_imgsrc_orig, delay_DS_imgsrc_orig, dir_DS_imgsrc_orig] = imgsrc_calculate(pos_list_orig,pos_src_orig,ord_imgsrc,dim_room,fs);
% view([0 90])
% ylim([0 dim_room(2)])
% xlim([0 dim_room(1)])
% pbaspect([dim_room(1) dim_room(2) 1])

%{
% debug - how well matched is the ISM and the measurement?

delay_diff_DS = delay_DS_orig_meas - delay_DS_imgsrc_orig;

delay_diff_threshold = 0.1; % threshold in ms before we need to correct it
delay_diff_threshold_samp = delay_diff_threshold/1000 * fs;

if abs(delay_diff_DS) > delay_diff_threshold_samp
    % if the measured and IS times are offset:
    disp(['Mismatch between time of direct sound of measurement and image source of ',num2str(delay_diff_DS/fs*1000,3),' ms. Counter-delaying imgsrc timings by ',num2str(round(delay_diff_DS)),' samples'])

    % counter-delay the image source times to nearest sample:
    delay_DS_imgsrc_orig = delay_DS_imgsrc_orig + round(delay_diff_DS);
    delay_imgsrc_orig = delay_imgsrc_orig + round(delay_diff_DS);
end
fig;
% subplot(3,1,1)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_orig(1:plotcutoff,1))))
hold on;
title(['Original. RMS = ',num2str(db(rms(srir_orig(:,1))),4),'dB'])
disp(['Original. RMS = ',num2str(db(rms(srir_orig(:,1))),4),'dB'])

text(round(delay_DS_imgsrc_orig)/fs,db(abs(srir_orig(round(delay_DS_imgsrc_orig),1))),'x DS','Color','red','Clipping','on')
for i = 1:size(delay_imgsrc_orig,2)
    text(round(delay_imgsrc_orig(1,i))/fs,db(abs(srir_orig(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)],'Color','red','Clipping','on')
end
ylim([-40 0])
% plot(plot_time_vector,db(abs(srir_tar(1:end,1))))
% xticks(1:500:plotcutoff)
% xticklabels = string(plot_time_vector(1:500:plotcutoff));
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])




%}

% target
[dist_imgsrc_tar, delay_imgsrc_tar, dir_imgsrc_tar, dist_DS_imgsrc_tar, delay_DS_imgsrc_tar, dir_DS_imgsrc_tar] = imgsrc_calculate(pos_list_tar,pos_src_tar,ord_imgsrc,dim_room,fs);
% view([0 90])
% ylim([0 dim_room(2)])
% xlim([0 dim_room(1)])
% pbaspect([dim_room(1) dim_room(2) 1])

% Sort arrivals based on lowest distances (mean of original and target distances)
% This ensures the first (most important) arrivals are processed first. 
[~,idx]=sort(mean([dist_imgsrc_orig; dist_imgsrc_tar]));

% If imgsrc 2nd order, it seems to duplicate the direct sound path. Remove the index of the shortest path to get rid:
if ord_imgsrc > 1 % if 2nd order or higher
    idx = idx(2:end);
end

% Sort
dist_imgsrc_orig = dist_imgsrc_orig(idx);
delay_imgsrc_orig = delay_imgsrc_orig(idx);
dir_imgsrc_orig = dir_imgsrc_orig(:,idx);

dist_imgsrc_tar = dist_imgsrc_tar(idx);
delay_imgsrc_tar = delay_imgsrc_tar(idx);
dir_imgsrc_tar = dir_imgsrc_tar(:,idx);


%% Compare measured and image source direct sound

% Get difference in direction of direct sound between measured and image
% source:
dir_diff_DS = rad2deg(angdiff(deg2rad(dir_DS_meas_orig), deg2rad(dir_DS_imgsrc_orig)));

dir_diff_threshold = 30; % threshold in degrees before we need to rotate it
% if the difference in direction of direct sound is greater than dir_diff_threshold:
if abs(dir_diff_DS(1,1))> dir_diff_threshold
    % round to nearest degrees that passes the dir_diff_threshold. Don't
    % want to be doing adjustments for only tiny differences!
    dir_diff_DS_rnd = dir_diff_threshold*round(dir_diff_DS/dir_diff_threshold);
    disp(['Mismatch between rotation of measurement and image source. Azimuth error = ',num2str(dir_diff_DS(1,1),3),'degrees. Counter-rotating by ',num2str(dir_diff_DS_rnd(1,1)),' degrees'])

    % rotate the image source directions:
    dir_DS_imgsrc_orig = dir_DS_imgsrc_orig-dir_diff_DS_rnd;
    dir_DS_imgsrc_tar = dir_DS_imgsrc_tar-dir_diff_DS_rnd;

    dir_imgsrc_orig = dir_imgsrc_orig-dir_diff_DS_rnd;
    dir_imgsrc_tar = dir_imgsrc_tar-dir_diff_DS_rnd;
end

%% now look at time differences between measured and image source direct sound:
delay_diff_DS = delay_DS_orig_meas - delay_DS_imgsrc_orig;

delay_diff_threshold = 0.1; % threshold in ms before we need to correct it
delay_diff_threshold_samp = delay_diff_threshold/1000 * fs;

if abs(delay_diff_DS) > delay_diff_threshold_samp
    % if the measured and IS times are offset:
    disp(['Mismatch between time of direct sound of measurement and image source of ',num2str(delay_diff_DS/fs*1000,3),' ms. Counter-delaying imgsrc timings by ',num2str(round(delay_diff_DS)),' samples'])

    % counter-delay the image source times to nearest sample:
    delay_DS_imgsrc_orig = delay_DS_imgsrc_orig + round(delay_diff_DS);
    delay_imgsrc_orig = delay_imgsrc_orig + round(delay_diff_DS);
    delay_DS_imgsrc_tar = delay_DS_imgsrc_tar + round(delay_diff_DS);
    delay_imgsrc_tar = delay_imgsrc_tar + round(delay_diff_DS);
end

%% Work out transform difference between original and target positions - to go from original to target
% for direct sound and reflections, and for direction, delay and loudness

% Direction / direction of arrival
% e.g. if orig -> tar means an arrival moves to right, need to move original to right
dir_diff_DS_transform = rad2deg(angdiff(deg2rad(dir_DS_imgsrc_orig), deg2rad(dir_DS_imgsrc_tar)));
dir_diff_transform = rad2deg(angdiff(deg2rad(dir_imgsrc_orig), deg2rad(dir_imgsrc_tar)));

% Delay / time of arrival
% e.g. if orig -> tar means arrival gets further away, delay should increase
delay_diff_DS_transform = round(delay_DS_imgsrc_tar - delay_DS_imgsrc_orig); % in samples
delay_diff_transform = round(delay_imgsrc_tar - delay_imgsrc_orig); % in samples

% Distance / gain - use the ratio of old to new distance to calculate gain
% e.g. if orig -> tar means an arrival increases in distance, level should decrease
dist_diff_DS_transform = (dist_DS_imgsrc_orig ./ dist_DS_imgsrc_tar);
dist_diff_transform = (dist_imgsrc_orig ./ dist_imgsrc_tar);

plotcutoff = fs/20;
%% Put x on the direct amplitude plot to see where our detected arrivals are

%{
if flag_subspace_decomp
    fig;hold on; plot(db(abs(srir_resid(1:plotcutoff,1)))); plot(db(abs(srir_direct(1:plotcutoff,1))));
    text(round(delay_DS_imgsrc_orig),db(abs(srir_direct(round(delay_DS_imgsrc_orig),1))),'x DS')
    for i = 1:size(delay_imgsrc_orig,2)
        text(round(delay_imgsrc_orig(1,i)),db(abs(srir_direct(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)],'Clipping','on')
    end
else
    fig;hold on; plot(db(abs(srir_orig(1:plotcutoff,1)))); %plot(db(abs(srir_direct(:,1))));
    text(round(delay_DS_imgsrc_orig),db(abs(srir_orig(round(delay_DS_imgsrc_orig),1))),'x DS')
    for i = 1:size(delay_imgsrc_orig,2)
        text(round(delay_imgsrc_orig(1,i)),db(abs(srir_orig(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)],'Clipping','on')
    end
end
xlim([0 plotcutoff])
ylim([-40 0])
%}
%% ================ EXTRAPOLATION ================
% Window out the arrivals from the SRIR, rotate + delay + gain, and add back in

gainpowerDS = 1; % 1 = linear, 2 = inverse square law;
gainpower = 1; % 1 = linear, 2 = inverse square law;

% create the new SRIR:
if flag_subspace_decomp
    srir_arrival_removed = srir_direct;
    srir_unAlt = srir_direct;
else
    srir_arrival_removed = srir_orig;
    srir_unAlt = srir_orig;
end

figure
% subplot(3,1,1)
plot(srir_unAlt(1:fs/10,1))
% hold on
% ylim([-1 1])

% % % % % % % % EARLY REFLECTIONS:
for i= 1:length(delay_imgsrc_orig) % for each image source

    % text(delay_imgsrc_orig(i)-delay_arrival_initial_samp,srir_arrival_removed(delay_imgsrc_orig(i)-delay_arrival_initial_samp,1),['X',num2str(i)],'Color','green')
    % text(delay_imgsrc_orig(i)-delay_arrival_initial_samp+arrival_dur_ms/1000*fs,srir_arrival_removed(delay_imgsrc_orig(i)-delay_arrival_initial_samp+arrival_dur_ms/1000*fs,1),['X',num2str(i)],'Color','red')

    % remove arrival
    % Two options here. Remove from the original (unaltered) SRIR, or remove
    % from the new SRIR. I think from unaltered is best. 
    % [srir_arrival,srir_arrival_removed] = srir_arrival_remove(srir_arrival_removed,srir_arrival_removed,fs,delay_imgsrc_orig(i)-delay_arrival_initial_samp,arrival_dur_ms);
    [srir_arrival,srir_arrival_removed] = srir_arrival_remove(srir_arrival_removed,srir_unAlt,fs,delay_imgsrc_orig(i)-delay_arrival_initial_samp,arrival_dur_ms);

    % rotate
    % Two options here. Both should in theory be the same, if the image source is
    % the same as the measured. Better results found getting doa of actual measurement rather than the IS.
    % % First option: based on doa of image source
    % srir_arrival_rot = rotateHOA_N3D(srir_arrival,(dir_diff_transform(1,i)),(dir_diff_transform(2,i)),0);
    % % Second option: based on doa of measurement
    dir_arr_orig = doa_iv(srir_arrival);
    dir_arr_orig = dir_imgsrc_orig(:,i);
    dir_arr_tar = dir_imgsrc_tar(:,i);
    dir_rot = rad2deg(angdiff(deg2rad(dir_arr_orig), deg2rad(dir_arr_tar)));
    srir_arrival_rot = rotateHOA_N3D(srir_arrival,dir_rot(1),dir_rot(2),0);

    % doa_iv(srir_arrival_rot)

    % dir_diff_transform(:,i)
    % gain and add back in at new time
    srir_arrival_removed(delay_imgsrc_tar(i)-delay_arrival_initial_samp+1:delay_imgsrc_tar(i)-delay_arrival_initial_samp+size(srir_arrival,1),:) ...
        = srir_arrival_removed(delay_imgsrc_tar(i)-delay_arrival_initial_samp+1:delay_imgsrc_tar(i)-delay_arrival_initial_samp+size(srir_arrival,1),:) + ...
        srir_arrival_rot*dist_diff_transform(i).^gainpower;

    %     figure(5)
    % plot(srir_arrival_removed(1:fs/20,1))
    % 
    % figure(11)
    % plot(srir_arrival_rot(:,1))
end


% % % % % % % % DIRECT SOUND:
[srir_arrival,srir_arrival_removed] = srir_arrival_remove(srir_arrival_removed,srir_unAlt,fs,delay_DS_imgsrc_orig-delay_arrival_initial_samp_DS,arrival_dur_ms_DS);
% doa_iv(srir_arrival)
% fig;plot(srir_arrival(:,1))
% dir_arr_orig = doa_iv(srir_arrival);
dir_arr_orig = dir_DS_imgsrc_orig;
dir_arr_tar = dir_DS_imgsrc_tar;
dir_rot = rad2deg(angdiff(deg2rad(dir_arr_orig), deg2rad(dir_arr_tar)));

% rotate it
% srir_arrival_rot = rotateHOA_N3D(srir_arrival,(dir_diff_DS_transform(1)),(dir_diff_DS_transform(2)),0);
% doa_iv(srir_arrival_rot)
srir_arrival_rot = rotateHOA_N3D(srir_arrival,dir_rot(1),dir_rot(2),0); % adhoc

% set to zero anything up to the new direct sound arrival time
% (important for when the target source is further than the original
% source!) 
srir_arrival_removed(1:delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+1,:) = 0;
srir_resid(1:delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+1,:) = 0;

% gain and add back in at new time
srir_arrival_removed(delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+1:delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+size(srir_arrival,1),:) ...
    = srir_arrival_removed(delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+1:delay_DS_imgsrc_tar-delay_arrival_initial_samp_DS+size(srir_arrival,1),:) + ...
    srir_arrival_rot*dist_diff_DS_transform.^gainpowerDS;

% figure(5)
% subplot(3,1,2)
% plot(srir_arrival_removed(1:fs/20,1))
% ylim([-1 1])

% figure(11)
% plot(srir_arrival_rot(:,1))


% figure(6)
% subplot(3,1,3)
hold on
plot(srir_arrival_removed(1:fs/10,1))
ylim([-1 1])
legend({'Original','Extrapolated'})

% combine with old residual for subspace decomp
if flag_subspace_decomp
    srir_new = srir_arrival_removed + srir_resid;
else
    srir_new = srir_arrival_removed;
end


%% ================ EVALUATION ================
%% Compare new and old SRIRs
% Time-domain:

fig;
subplot(3,1,1)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_orig(1:plotcutoff,1))))
hold on;
title(['Original. RMS = ',num2str(db(rms(srir_orig(:,1))),4),'dB'])
disp(['Original. RMS = ',num2str(db(rms(srir_orig(:,1))),4),'dB'])

text(round(delay_DS_imgsrc_orig)/fs,db(abs(srir_orig(round(delay_DS_imgsrc_orig),1))),'x DS','Color','red','Clipping','on')
for i = 1:size(delay_imgsrc_orig,2)
    text(round(delay_imgsrc_orig(1,i))/fs,db(abs(srir_orig(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)],'Color','red','Clipping','on')
end
ylim([-40 0])
% plot(plot_time_vector,db(abs(srir_tar(1:end,1))))
% xticks(1:500:plotcutoff)
% xticklabels = string(plot_time_vector(1:500:plotcutoff));
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])

subplot(3,1,2)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_new(1:plotcutoff,1))))
hold on
title(['Extrapolated. RMS = ',num2str(db(rms(srir_new(:,1))),4),'dB'])
disp(['Extrapolated. RMS = ',num2str(db(rms(srir_new(:,1))),4),'dB'])
% plot(db(abs(srir_new(1:fs/10,1))))
text(round(delay_DS_imgsrc_tar)/fs,db(abs(srir_new(round(delay_DS_imgsrc_tar),1))),'+ DS','Color','green','Clipping','on')
for i = 1:size(delay_imgsrc_orig,2)
    text(round(delay_imgsrc_tar(1,i))/fs,db(abs(srir_new(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)],'Color','green','Clipping','on')
end
text(round(delay_DS_imgsrc_orig)/fs,db(abs(srir_new(round(delay_DS_imgsrc_orig),1))),'- DS','Color','red','Clipping','on')
for i = 1:size(delay_imgsrc_orig,2)
    text(round(delay_imgsrc_orig(1,i))/fs,db(abs(srir_new(round(delay_imgsrc_orig(1,i)),1))),['- ', num2str(i)],'Color','red','Clipping','on')
end
ylim([-40 0])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])

subplot(3,1,3)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_tar(1:plotcutoff,1))))
hold on;
title(['Target. RMS = ',num2str(db(rms(srir_tar(:,1))),4),'dB'])
disp(['Target. RMS = ',num2str(db(rms(srir_tar(:,1))),4),'dB'])
text(round(delay_DS_imgsrc_tar)/fs,db(abs(srir_tar(round(delay_DS_imgsrc_tar),1))),'+ DS','Color','green','Clipping','on')
for i = 1:size(delay_imgsrc_orig,2)
    text(round(delay_imgsrc_tar(1,i))/fs,db(abs(srir_tar(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)],'Color','green','Clipping','on')
end
ylim([-40 0])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])

%% SIMPLER VERSION:
%{
fig;
% subplot(3,1,1)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_orig(1:plotcutoff,1))))
hold on;
% title(['Original. RMS = ',num2str(db(rms(srir_orig(:,1))),4),'dB'])
% text(round(delay_DS_imgsrc_orig)/fs,db(abs(srir_orig(round(delay_DS_imgsrc_orig),1))),'x DS','Color','red','Clipping','on')
% for i = 1:size(delay_imgsrc_orig,2)
%     text(round(delay_imgsrc_orig(1,i))/fs,db(abs(srir_orig(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)],'Color','red','Clipping','on')
% end
ylim([-40 0])
% plot(plot_time_vector,db(abs(srir_tar(1:end,1))))
% xticks(1:500:plotcutoff)
% xticklabels = string(plot_time_vector(1:500:plotcutoff));
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])
pbaspect([3 1 1]);
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_orig.pdf'); end

fig
% subplot(3,1,2)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_new(1:plotcutoff,1))))
hold on
% title(['Extrapolated. RMS = ',num2str(db(rms(srir_new(:,1))),4),'dB'])
% % plot(db(abs(srir_new(1:fs/10,1))))
% text(round(delay_DS_imgsrc_tar)/fs,db(abs(srir_new(round(delay_DS_imgsrc_tar),1))),'+ DS','Color','green','Clipping','on')
% for i = 1:size(delay_imgsrc_orig,2)
%     text(round(delay_imgsrc_tar(1,i))/fs,db(abs(srir_new(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)],'Color','green','Clipping','on')
% end
% text(round(delay_DS_imgsrc_orig)/fs,db(abs(srir_new(round(delay_DS_imgsrc_orig),1))),'- DS','Color','red','Clipping','on')
% for i = 1:size(delay_imgsrc_orig,2)
%     text(round(delay_imgsrc_orig(1,i))/fs,db(abs(srir_new(round(delay_imgsrc_orig(1,i)),1))),['- ', num2str(i)],'Color','red','Clipping','on')
% end
ylim([-40 0])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])
pbaspect([3 1 1]);
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_extrap.pdf');end

fig
% subplot(3,1,3)
plot(plot_time_vector(1:plotcutoff),db(abs(srir_tar(1:plotcutoff,1))))
hold on;
% title(['Target. RMS = ',num2str(db(rms(srir_tar(:,1))),4),'dB'])
% text(round(delay_DS_imgsrc_tar)/fs,db(abs(srir_tar(round(delay_DS_imgsrc_tar),1))),'+ DS','Color','green','Clipping','on')
% for i = 1:size(delay_imgsrc_orig,2)
%     text(round(delay_imgsrc_tar(1,i))/fs,db(abs(srir_tar(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)],'Color','green','Clipping','on')
% end
ylim([-40 0])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])
pbaspect([3 1 1]);
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_target.pdf');end


% % this was with direct and residual plotted separately
% subplot(3,1,1)
% title('Original')
% plot(db(abs(srir_resid(1:plotcutoff,1))))
% ylim([-50 0])
% hold on;
% plot(db(abs(srir_direct(1:plotcutoff,1))))
%
% text(round(delay_DS_imgsrc_orig),db(abs(srir_direct(round(delay_DS_imgsrc_orig),1))),'x DS')
% for i = 1:size(delay_imgsrc_orig,2)
% text(round(delay_imgsrc_orig(1,i)),db(abs(srir_direct(round(delay_imgsrc_orig(1,i)),1))),['x ', num2str(i)])
% end
%
% subplot(3,1,2)
% title('Extrap')
% plot(db(abs(srir_resid(1:plotcutoff,1))))
% hold on
% plot(db(abs(srir_arrival_removed(1:plotcutoff,1))))
% % plot(db(abs(srir_new(1:fs/10,1))))
% text(round(delay_DS_imgsrc_tar),db(abs(srir_arrival_removed(round(delay_DS_imgsrc_tar),1))),'+ DS')
% for i = 1:size(delay_imgsrc_orig,2)
% text(round(delay_imgsrc_tar(1,i)),db(abs(srir_arrival_removed(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)])
% end
%
% text(round(delay_DS_imgsrc_orig),db(abs(srir_arrival_removed(round(delay_DS_imgsrc_orig),1))),'- DS')
% for i = 1:size(delay_imgsrc_orig,2)
% text(round(delay_imgsrc_orig(1,i)),db(abs(srir_arrival_removed(round(delay_imgsrc_orig(1,i)),1))),['- ', num2str(i)])
% end
%
% ylim([-50 0])
%
% subplot(3,1,3)
% title('Target')
% plot(db(abs(srir_tar(1:plotcutoff,1))))
% ylim([-50 0])
% hold on;
% text(round(delay_DS_imgsrc_tar),db(abs(srir_tar(round(delay_DS_imgsrc_tar),1))),'+ DS')
% for i = 1:size(delay_imgsrc_orig,2)
% text(round(delay_imgsrc_tar(1,i)),db(abs(srir_tar(round(delay_imgsrc_tar(1,i)),1))),['+ ', num2str(i)])
% end

%if save_figs;  exportgraphics(gcf, 'srir_extrap_TD.pdf'); end

%}

%% Compare extrapolated SRIR to target SRIR
% %{
% in same plot - see what's moved
fig
% subplot(2,1,1)
hold on
ylim([-30 0])
% plot(db(abs(srir_direct(1:plotcutoff,1))))
% plot(db(abs(srir_arrival_removed(1:plotcutoff,1))))
% plot(db(abs(srir_tar(1:plotcutoff,1))))
plot(plot_time_vector(1:plotcutoff),db(abs(srir_orig(1:plotcutoff,1))))
plot(plot_time_vector(1:plotcutoff),db(abs(srir_new(1:plotcutoff,1))))
% plot(db(abs(srir_tar(1:plotcutoff,1))))
legend({'Original','Extrapolated','Target'})
% xlim([1 plotcutoff])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])
pbaspect([3 1 1]);
box on
grid on
% title(['Original']);
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_comparison_orig.pdf');end

fig
% subplot(2,1,2)
hold on
ylim([-30 0])
% plot(db(abs(srir_orig(1:plotcutoff,1))))
plot(plot_time_vector(1:plotcutoff),db(abs(srir_tar(1:plotcutoff,1))))
plot(plot_time_vector(1:plotcutoff),db(abs(srir_new(1:plotcutoff,1))))
% plot(db(abs(srir_tar(1:plotcutoff,1))))
legend({'Target','Extrapolated'})
% xlim([1 plotcutoff])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
xlim([0 max(plot_time_vector(1:plotcutoff))])
pbaspect([3 1 1]);
box on
grid on
if save_figs; exportgraphics(gcf, 'srir_extrap_TD_comparison_extrap.pdf');end
% exportgraphics(gcf, 'srir_extrap_TD_comparison.pdf');
%}

%% Binaural render
    %render and play IR_original(9) IR_generate IR_record
    % 1.S3L3 2.generete 3.target recoeding
    % doBinRenderIR=1;
    % if doBinRenderIR
        % IR_Record=squeeze(sofa1.Data.IR(idx_LS_srcrec_orig,:,:)) ;
        % srir_orig=squeeze(sofa1.Data.IR(9,:,:));
        brir_orig=binSound(srir_orig,SH_ambisonic_binaural_decoder);
        brir_new=binSound(srir_new,SH_ambisonic_binaural_decoder);
        brir_tar=binSound(srir_tar,SH_ambisonic_binaural_decoder);
        % gaptime=zeros(0.3*fs,2);
        %soundsc([binIR_original;gaptime;binIR_generate;gaptime;binIR_record],Fs);
        %soundsc([binIR_generate;gaptime;binIR_record],Fs);% only last 2 sound

%% Calculate colouration

brir_o = reshape(brir_orig,[],1,2);
brir_n = reshape(brir_new,[],1,2);
brir_t = reshape(brir_tar,[],1,2);

settings.smGL1 = 0;

pbc2_orig = mckenzie2025(brir_o,brir_t,settings);
pbc2_new = mckenzie2025(brir_n,brir_t,settings);
% pbc2_tar = mckenzie2025(brir_tar,brir_tar);

% Predicted binaural colouration

disp(['PBC (orig -> tar) = ',num2str(pbc2_orig,3)]);
disp(['PBC (new -> tar) = ',num2str(pbc2_new,3)]);


%% Frequency plot (W channel)
%{
fig
freqplot_smooth(srir_orig(:,1),fs,2)
hold on
freqplot_smooth(srir_new(:,1),fs,2)
freqplot_smooth(srir_tar(:,1),fs,2)
% legend({['Original. ','PBC = ',num2str(pbc2_orig,3)],['Extrapolated. ','PBC = ',num2str(pbc2_new,3)],'Target'})
legend({['Original'],['Extrapolated'],'Target'})
pbaspect([2 1 1]);
ylim([-5 25])
xlim([40 20000])

%}
% exportgraphics(gcf, 'srir_extrap_FD_comparison.pdf');

%% Frequency plot (W channel) comparison
%{
fig
freqplot_smooth(srir_orig(:,1),fs)
hold on
freqplot_smooth(srir_new(:,1),fs)
% freqplot_smooth(srir_tar(1:plotcutoff,1),fs)
% legend({['Original. ','PBC = ',num2str(pbc2_orig,3)],['Extrapolated. ','PBC = ',num2str(pbc2_new,3)],'Target'})
legend({['Original'],['Extrapolated'],'Target'})
% pbaspect([2 1 1]);
pbaspect([3 1 1]);

ylim([-4 23])
xlim([40 20000])

if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_orig.pdf');end

fig
freqplot_smooth(srir_tar(:,1),fs)
hold on
freqplot_smooth(srir_new(:,1),fs)
% legend({['Original. ','PBC = ',num2str(pbc2_orig,3)],['Extrapolated. ','PBC = ',num2str(pbc2_new,3)],'Target'})
legend({'Target','Extrapolated'})
% pbaspect([2 1 1]);
pbaspect([3 1 1]);

ylim([-4 23])
xlim([40 20000])
if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_extrap.pdf');end

%}
%% Frequency plot (binaural) comparison

fig
subplot(1,2,1)
freqplot_smooth(brir_orig(:,1),fs)
hold on
freqplot_smooth(brir_new(:,1),fs)
pbaspect([1.5 1 1]);
ylim([0 35])
xlim([40 20000])
title('Left')

subplot(1,2,2)
freqplot_smooth(brir_orig(:,2),fs)
% freqplot_smooth(brir_new(:,1),fs)
hold on
freqplot_smooth(brir_new(:,2),fs)
% legend({['Original'],['Extrapolated'],'Target'})
ylabel('')
yticklabels('')
% legend({'Original (left)','Original (right)','Extrapolated (left)','Extrapolated (right)'})
pbaspect([1.5 1 1]);
ylim([0 35])
xlim([40 20000])
legend({['Original'],['Extrapolated'],'Target'},'location','southwest')
title('Right')

aa=subplot(122);
aa.Position(1)=0.5;

% fig
% freqplot_smooth(brir_tar(:,1),fs)
% hold on
% freqplot_smooth(brir_tar(:,2),fs)
% 
% freqplot_smooth(brir_new(:,1),fs)
% freqplot_smooth(brir_new(:,2),fs)
% 
% legend({'Target (left)','Target (right)','Extrapolated (left)','Extrapolated (right)'})
% pbaspect([3 1 1]);
% 
% ylim([0 35])
% xlim([40 20000])

if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_orig_bin.pdf');end

fig
subplot(1,2,1)
freqplot_smooth(brir_tar(:,1),fs)
hold on
freqplot_smooth(brir_new(:,1),fs)
pbaspect([1.5 1 1]);
ylim([0 35])
xlim([40 20000])
title('Left')

subplot(1,2,2)
freqplot_smooth(brir_tar(:,2),fs)
% freqplot_smooth(brir_new(:,1),fs)
hold on
freqplot_smooth(brir_new(:,2),fs)
% legend({['Original'],['Extrapolated'],'Target'})
ylabel('')
yticklabels('')
% legend({'Original (left)','Original (right)','Extrapolated (left)','Extrapolated (right)'})
pbaspect([1.5 1 1]);
ylim([0 35])
xlim([40 20000])
legend({'Target','Extrapolated'},'location','southwest')
title('Right')

aa=subplot(122);
aa.Position(1)=0.5;


if save_figs; exportgraphics(gcf, 'srir_extrap_FD_comparison_tar_bin.pdf');end



%% Plot horizontal DoA

srir_combined(:,:,1) =  srir_orig;
srir_combined(:,:,2) =  srir_new;
srir_combined(:,:,3) =  srir_tar;

% subplot(3,1,1);title('Original')
fig;
[doa_hor,doa_hor_p,P_pwd] = plot_doa_horiz(srir_combined,fs);
% xlim([0.5 size(srir_combined,3)+0.5]);
% xticks(1:3);
yticklabels({'Original','Extrapolated','Target'})
% subplot(3,1,2);title('Extrapolated')
% [doa_extrap,doa_extrap_p] = plot_doa_horiz(srir_new,fs);
% subplot(3,1,3);title('Target')
% [doa_tar,doa_tar_p] = plot_doa_horiz(srir_tar,fs);

set(gca,'FontSize',12)

pbaspect([3 1 1]);
if save_figs; exportgraphics(gcf, 'srir_extrap_hor_doa.pdf');end

% grid minor

% pause;


%% Plot DoA 2D 
%{
% normalized_doa_est_P = P_pwd./ max(P_pwd,[],1);
% normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;
% 
% normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;
% 
% P_pwd_n = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);
plot_thresh = -15;
normalized_doa_est_P = P_pwd./ max(P_pwd,[],1);
normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;

normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;

% P_pwd_n = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);
P_pwd_n= normalized_doa_est_P_dB;


degreeResolution = 2;

% P_pwd_n = doa_hor_p;

h=figure;
x = 1:length(srir_combined(1,1,:));
y = -180:degreeResolution:180-degreeResolution;

% xlin = linspace(min(x),max(x),180*2);
% ylin = linspace(min(y),max(y),90*2);
% [X,Y] = meshgrid(xlin,ylin);
[X,Y] = meshgrid(x,y);

Z = griddata(x,y,P_pwd_n(:,:),X,Y,'cubic');



surf(Y,X,Z,'EdgeColor','none')
axis tight;
xlabel('Azimuth (°)');
% ylabel('RIR number');
% 
yticklabels({'Original','Extrapolated','Target'})

shading interp
view ([0 90])
% % s.EdgeColor = 'none';
% 
set(gca, 'ActivePositionProperty' , 'position');

originalSize = get(gca, 'Position');
originalOuterBoundary = get(gca, 'outerposition');

set(gca, 'Position', originalSize);
set(gca, 'outerposition', originalOuterBoundary);
% 
colormap(flipud(bone))
set(gca, 'XTick', -150:75:150); %%%% should that be -180:75:180   ???
xlim([-180 180]);
% ylim([-90 90]);

% m = pbaspect('mode');
% m = 'manual';
pbaspect([4 1 1])
% daspect([1 0.6 1]);
set(gca, 'fontsize', 12);
set(gcf, 'Color', 'w');
box on



% view([90 -90])
c2 = colorbar;
c2.Label.String = 'Normalised Power';
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')

%}

%% Plot DoA 3D
figure;

% subplot(3,1,1);
[P_pwd,~,~,grid_dirs] = get_pwd(srir_orig,fs);
heatmap_plot(rad2deg(grid_dirs(:,1)),rad2deg(grid_dirs(:,2)),P_pwd);
colormap(flipud(bone))
% clim([0 2])
clim([0 3])

% clim([min(P_pwd) max(P_pwd)])
% title('Original')
set(gca,'FontSize',11)
k = colorbar;
xlabel(k,'Normalised power (dB)');
%caxis([plot_thresh 0]);
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_original.pdf');end

% max(P_pwd)

fig
% subplot(3,1,2);
[P_pwd,~,~,grid_dirs] = get_pwd(srir_new,fs);
heatmap_plot(rad2deg(grid_dirs(:,1)),rad2deg(grid_dirs(:,2)),P_pwd);
colormap(flipud(bone))
% clim([0 2])
clim([0 3])

% clim([min(P_pwd) max(P_pwd)])

% title('Extrapolated')
set(gca,'FontSize',11)
% box on

k = colorbar;
xlabel(k,'Normalised power (dB)');
%caxis([plot_thresh 0]);
% min(P_pwd)
% max(P_pwd)
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_extrap.pdf');end

fig
% subplot(3,1,3);
[P_pwd,~,~,grid_dirs] = get_pwd(srir_tar,fs);
heatmap_plot(rad2deg(grid_dirs(:,1)),rad2deg(grid_dirs(:,2)),P_pwd);
colormap(flipud(bone))
clim([0 3])
% clim([min(P_pwd) max(P_pwd)])

% title('Target')
set(gca,'FontSize',11)
k = colorbar;
xlabel(k,'Normalised power (dB)');
%caxis([plot_thresh 0]);
% min(P_pwd)
% max(P_pwd)
if save_figs; exportgraphics(gcf, 'srir_extrap_doa_target.pdf');end

%% Very simple time-domain plot
%{
fig
subplot(3,1,1)
plot_time_vector = 0:1/fs:(size(srir_orig,1)/fs-1/fs);
plot(plot_time_vector,db(abs(srir_orig(1:end,1))))
ylim([-60 0])
ylabel('Magnitude (dB)')
xlabel('Time (s)')
title('Original')

subplot(3,1,2)
plot(plot_time_vector,db(abs(srir_new(1:end,1))))
ylim([-60 0])
ylabel('Magnitude (dB)')
title('Extrapolated')
xlabel('Time (s)')

subplot(3,1,3)
plot(plot_time_vector,db(abs(srir_tar(1:end,1))))
ylim([-60 0])
ylabel('Magnitude (dB)')
title('Target')
xlabel('Time (s)')



%}


% legend({'Original','Extrapolated','Target'})

% doa = doa_iv(srir_orig(1:450,:))
% doa = doa_iv(srir_new(1:450,:))
% doa = doa_iv(srir_tar(1:450,:))

%% LISTEN

soundsc([srir_orig(:,1);srir_new(:,1);srir_tar(:,1)],fs)

%%
% SourcerPoint_generate=sofa1.SourcePosition(idx_LS_srcrec_orig,:);
% ListenerPoint_generate=sofa1.ListenerPosition(idx_LS_srcrec_orig,:);% the controlRecordNum recording it target
% arrival_time_original=delay_DS_orig_meas(9);
% %Am_original=directSoundValure(9);
% % cut
% directSound=srir_orig(:,:,arrival_time_original-directSoundCutLeft:arrival_time_original+directSoundCutRight);
%
% %fit distance - amiptude
% %use rms
% Nface=15;%only facing 15
% rmsDirectSoundValure=zeros(Nface,1);
%
% for n=1:Nface
%     readIR_forFit=sofa1.Data.IR(n,:,:);
%     arrival_time_forFit=firstReflectionIndex(n);
%     ds=readIR_forFit(:,1,arrival_time_forFit-directSoundCutLeft:arrival_time_forFit+directSoundCutRight);
%     rmsDirectSoundValure(n)=rms(squeeze(ds));
% end
% fun = @(x)sseval(x,dist_srcrec(1:Nface,1),rmsDirectSoundValure);
% x0 = rand(2,1);
% bestx = fminsearch(fun,x0);
%
% pic2=1;
% if pic2
%     figure;% check the fit
%     scatter(dist_srcrec,rmsDirectSoundValure);
%     for n=1:N
%         text(dist_srcrec(n),rmsDirectSoundValure(n),num2str(n));
%     end
%     hold on;
%     dislist=0:0.2:10;
%     plot(dislist,bestx(1)*exp(-(bestx(2))*dislist));
%     xlabel('distace/m')
%     ylabel('rms Amplitude of Direct Sound')
%
% end
%
% %compute the time and Amiplitude
% distance_generate=sqrt(sum((SourcerPoint_generate-ListenerPoint_generate).^2));
% distance_origianl=sqrt(sum((pos_src_orig-pos_list_orig).^2));
% arrival_time_generate=floor(distance_generate/speed);       %  in samples
% Am_generate=bestx(1)*exp(-(bestx(2))*distance_generate);
% Am_original=bestx(1)*exp(-(bestx(2))*distance_origianl);
% directSound=directSound*(Am_generate/Am_original); %scale
%
% % re-join
% A=zeros(25,arrival_time_generate-directSoundCutLeft-1);
% B=squeeze(directSound);         %direct sound
% C=squeeze(srir_orig(:,:,arrival_time_generate+directSoundCutRight+1:arrival_time_generate+directSoundCutRight+earlyRefCutLength));%early ref
% %C=C*(Am_generate/Am_original)/0.25;%scale C
% D=squeeze(srir_orig(:,:,arrival_time_generate+directSoundCutRight+earlyRefCutLength+1:end));%reverb
%
% IR_scale=[A,B,C,D];
%
% %% plot and listen
% %plot after scale
% pic3=0;
% if pic3
%     %sound(IR_scale(1,:),fs);
%     figure;
%
%     subplot(3,1,1);
%     plot(abs(squeeze(srir_orig(:,1,:))));
%     legend('IR original');
%     xlim([0 3000]);
%     ylim([0 0.5]);
%
%     subplot(3,1,2);
%     plot(abs(squeeze(IR_scale(1,:))));
%     legend('IR generate');
%     xlim([0 3000]);
%     ylim([0 0.5]);
%
%     subplot(3,1,3);
%     plot(abs(squeeze(sofa1.Data.IR(idx_LS_srcrec_orig,1,:))));
%     legend('IR record');
%     xlim([0 3000]);
%     ylim([0 0.5]);
% end
%
%
%
% %% rotate
%
% % compute the sph
% xyz1=pos_src_orig - pos_list_orig;
% xyz2=SourcerPoint_generate-ListenerPoint_generate;
% [az1,el1,r1]=cart2sph(xyz1(1),xyz1(2),xyz1(3));
% [az2,el2,r2]=cart2sph(xyz2(1),xyz2(2),xyz2(3));
% yaw=(az2-az1)/pi*180;
% pitch=(el2-el1)/pi*180;
% roll=0;
%
% %reshape hoasig
% hoasig=B.'; %B- reshape direct sound
%
% % rotate
% hoasig_rot = rotateHOA_N3D(hoasig, yaw, pitch, roll);
% B_rot=hoasig_rot.';
%
% %re-join
% IR_rot=[A,B_rot,C,D];
%
%
% %% plot and listen
% %plot after rotate
% pic3=1;
% if pic3
%     %sound(IR_scale(1,:),fs);
%     figure;
%
%     subplot(3,1,1);
%     plot(abs(squeeze(srir_orig(:,1,:))));
%     legend(['IR original ',num2str(3),num2str(3)]);
%     xlim([0 8000]);
%     ylim([0 0.5]);
%
%     subplot(3,1,2);
%     plot(abs(squeeze(IR_rot(1,:))));
%     legend('IR rot');
%     xlim([0 8000]);
%     ylim([0 0.5]);
%
%     subplot(3,1,3);
%     plot(abs(squeeze(sofa1.Data.IR(idx_LS_srcrec_orig,1,:))));
%     legend(['IR record ',num2str(idx_LS_src_orig),num2str(idx_LS_rec_orig)]);
%     xlim([0 8000]);
%     ylim([0 0.5]);
% end
%


%% render binaural and play

%run this to load decoder
if flag_render_binaural
    load_ambisonic_configuration


    % play all record
    playAllRecord=0;
    if playAllRecord
        for i=1:N/3
            speakerNum=3;
            recordNum=i*3-3+speakerNum;
            IR_test=squeeze(sofa1.Data.IR(recordNum,:,:)) ;
            waitingTime=zeros((i-1)*124223,2);
            binSound(IR_test,SH_ambisonic_binaural_decoder,Fs,waitingTime);
        end
    end

    %render and play IR_original(9) IR_generate IR_record
    % 1.S3L3 2.generete 3.target recoeding
    doBinRenderIR=1;
    if doBinRenderIR
        IR_Record=squeeze(sofa1.Data.IR(idx_LS_srcrec_orig,:,:)) ;
        srir_orig=squeeze(sofa1.Data.IR(9,:,:));

        binIR_original=binSound(srir_orig,SH_ambisonic_binaural_decoder);
        binIR_generate=binSound(IR_rot,SH_ambisonic_binaural_decoder);
        binIR_record=binSound(IR_Record,SH_ambisonic_binaural_decoder);
        gaptime=zeros(0.3*fs,2);
        %soundsc([binIR_original;gaptime;binIR_generate;gaptime;binIR_record],Fs);
        %soundsc([binIR_generate;gaptime;binIR_record],Fs);% only last 2 sound

    end

    %play controlRecord then generate IR conv with dry guitar
    % 1.S3L3 2.generete 3.target recoeding
    doBinRenderSONG=1;
    if doBinRenderSONG
        song_Dry=audioread('speechdirectsound_48.wav');
        [~,ch]=size(song_Dry);
        if ch==2
            song_Dry=(song_Dry(:,1)+song_Dry(:,2))/2;
        end

        clear binSong_recordIR;
        clear binSong_originalIR;
        clear binSong_generateIR;

        for i=1:2 % L and R
            binSong_originalIR(:,i)=conv(song_Dry,binIR_original(:,i));
            binSong_generateIR(:,i)=conv(song_Dry,binIR_generate(:,i));
            binSong_recordIR(:,i)=conv(song_Dry,binIR_record(:,i));
        end


        gaptime=zeros(0.5*fs,2);% gap 0.5s*fs
        soundsc([binSong_originalIR;gaptime;binSong_generateIR;gaptime;binSong_recordIR],fs);

        if ~exist('Max')
            Max=1;
        end


        audiowrite(['S',num2str(idx_LS_src_orig),'L',num2str(idx_LS_rec_orig),'_Original33','.wav'],binSong_originalIR/Max,fs);
        audiowrite(['S',num2str(idx_LS_src_orig),'L',num2str(idx_LS_rec_orig),'_Rotate','.wav'],binSong_generateIR/Max,fs);
        %audiowrite(['S',num2str(useSpeaker),'L',num2str(useListener),'_Record','.wav'],binSong_recordIR/Max,fs);

    end

end
%% declare functions

%Find  direct sound loc and val
function [locD,ValD,pks,lcs]=findDirectSound(ir)

highPassFilterFreq = 5000;
fs=48000;
[~,filtHi,~] = ambisonic_crossover(highPassFilterFreq,fs);
    ir = filter(filtHi,1,ir); % high pass filter
ir = circshift(ir,-floor(length(filtHi)/2));

absir=abs(ir);
absir = absir ./ max(absir);
noiseValuse=max(absir);
[pks,lcs]=findpeaks(absir,"MinPeakDistance",50,MinPeakHeight=0.05*noiseValuse);  %denosing,find peak
% Find the index of the peak representing the direct sound
[pk_max,ix_max] = max(pks);
locD=lcs(ix_max);
ValD=pks(ix_max);
end

% find A0 and alpha for distance - amipitude fit

function sse = sseval(x,tdata,ydata)
A = x(1);
lambda = x(2);
sse = sum((ydata - A*exp(-lambda*tdata)).^2);
end


function binIR=binSound(IR_test,SH_amb_bin_dec)% render to a bin audio and play it
ambisonic_soundscape=IR_test.';

binaural_ambisonic_render = zeros(length(SH_amb_bin_dec(1,:,1))+length(ambisonic_soundscape)-1,length(SH_amb_bin_dec(1,1,:)));

% convolve each channel of the encoded signal with the decoder signal and sum the result
for i = 1:length(SH_amb_bin_dec(:,1,1))

    binaural_ambisonic_render(:,1) = binaural_ambisonic_render(:,1) +  conv(SH_amb_bin_dec(i,:,1),ambisonic_soundscape(i,:))';
    binaural_ambisonic_render(:,2) = binaural_ambisonic_render(:,2) +  conv(SH_amb_bin_dec(i,:,2),ambisonic_soundscape(i,:) )';
end

% compare the two binaural decoders by listening to both binaural renders consecutively
binIR=binaural_ambisonic_render;
end


function [dist_imgsrc, delay_imgsrc, dir_imgsrc, dist_DS_imgsrc, delay_DS_imgsrc, dir_DS_imgsrc] = imgsrc_calculate(coord_rec,coord_src,ord_imgsrc,dim_room,fs)
% c = 329; % Speed of sound (m/s)
c = 343; % Speed of sound (m/s)
% c = 363; % Speed of sound (m/s)

h = figure;
% figure
plotRoom(dim_room,coord_rec,coord_src,h)

hold on

Lx = dim_room(1);
Ly = dim_room(2);
Lz = dim_room(3);

x = coord_src(1);
y = coord_src(2);
z = coord_src(3);


% JUST FIRST ORDER REFS:
if ord_imgsrc == 1
    xyz_src = [x -y -z;... % true 1st order
        -x  y -z;...
        -x  -y  z;...
        x -y -z;...
        -x  y -z;...
        -x  -y  z].';

    xyz_room =  [0 0 0;...
        0 0 0;...
        0 0 0;...
        2*Lx 0 0;...
        0 2*Ly 0;...
        0 0 2*Lz].';

    coords_imgsrc = xyz_room-xyz_src;

    for kk=1:size(xyz_src,2)
        coord_imgsrc=coords_imgsrc(:,kk);
        plot3(coord_imgsrc(1),coord_imgsrc(2),coord_imgsrc(3),"g*")
    end

else
    % Second or higher order IMG SRC
    xyz_src = [-x -y -z;... % matlab version
        -x -y  z;...
        -x  y -z;...
        -x  y  z;...
        x -y -z;...
        x -y  z;...
        x  y -z;...
        x  y  z].';
    % Increase the range to plot more images
    nVect = -(ord_imgsrc-1):(ord_imgsrc-1);
    lVect = -(ord_imgsrc-1):(ord_imgsrc-1);
    mVect = -(ord_imgsrc-1):(ord_imgsrc-1);

    coords_imgsrc_tot=[];

    for n = nVect
        for l = lVect
            for m = mVect
                xyz_diff = [n*2*Lx; l*2*Ly; m*2*Lz];
                coords_imgsrc = xyz_diff - xyz_src;
                coords_imgsrc_tot = [coords_imgsrc_tot  coords_imgsrc];
                for kk=1:size(xyz_src,2)
                    coord_imgsrc=coords_imgsrc(:,kk);
                    plot3(coord_imgsrc(1),coord_imgsrc(2),coord_imgsrc(3),"g*")
                    text(coord_imgsrc(1),coord_imgsrc(2),coord_imgsrc(3),[num2str(n),',',num2str(l),',',num2str(m),',',num2str(kk)])
                    plot3([coord_imgsrc(1);coord_rec(1)],[coord_imgsrc(2);coord_rec(2)],[coord_imgsrc(3);coord_rec(3)])
                end
            end
        end
    end
    coords_imgsrc = coords_imgsrc_tot;
end
% xlim([-3*Lx,3*Lx]);
% ylim([-3*Ly,3*Ly]);
% zlim([-3*Lz,3*Lz]);

% view([0 90])
%% Compute reflection delays and power

% n = 1;
% l = 1;
% m = 1;
% p = 1;



% Get the coordinates of the image.
% isourceCoord = [n*2*Lx; l*2*Ly; m*2*Lz] - sourceXYZ(:,p);


for i = 1:size(coords_imgsrc,2)
    coord_imgsrc = coords_imgsrc(:,i);

    % Calculate the delay (in samples) at which the contribution occurs.
    dist_imgsrc(i) = norm((coord_imgsrc(:)-coord_rec(:)),2); % distance from img src to rec
    delay_imgsrc(i) = round((fs/c).*dist_imgsrc(i)); % in samples from beginning, time = 0

    xyz_diff = coord_imgsrc-coord_rec';
    hyp = sqrt(xyz_diff(1)^2+xyz_diff(2)^2);
    elevation(i) = atan(xyz_diff(3)./(hyp+eps));
    azimuth(i) = atan2(xyz_diff(2),xyz_diff(1));

    % ImagePower = BX1.^abs(n-q(p)).*BY1.^abs(l-j(p)).*BZ1.^abs(m-k(p)).*BX2.^abs(n).*(BY2.^abs(l)).*(BZ2.^abs(m));
end

dir_imgsrc = [azimuth; elevation]*180/pi;

%%% tm might need converting ie plus 90 degrees azimuth?

% Also calculate the delay and direction of the direct sound
dist_DS_imgsrc = norm((coord_src(:)-coord_rec(:)),2); % distance from actual src to rec
delay_DS_imgsrc = round((fs/c).*dist_DS_imgsrc); % in samples from beginning, time = 0

xyz_diff = coord_src'-coord_rec';
hyp = sqrt(xyz_diff(1)^2+xyz_diff(2)^2);
elevation = atan(xyz_diff(3)./(hyp+eps));
azimuth = atan2(xyz_diff(2),xyz_diff(1));
dir_DS_imgsrc = [azimuth; elevation]*180/pi;





end


function [srir_arrival,srir_arrival_removed] = srir_arrival_remove(srir_in,srir_unAlt,fs,time_arrival_samp,arrival_dur_ms)

% arrival_dur_ms = 5; % how long to extract in ms
arrival_window_samp = [5 10]; % ramp on and ramp off - avoid clicks

arrival_window = ones(arrival_dur_ms/1000*fs,size(srir_in,2));
arrival_window(1:arrival_window_samp(1),:) = repmat(linspace(0,1,arrival_window_samp(1))',1,size(srir_in,2));
arrival_window(end-arrival_window_samp(2)+1:end,:) = repmat(linspace(1,0,arrival_window_samp(2))',1,size(srir_in,2));

srir_extract_window = srir_in * 0;
srir_extract_window(time_arrival_samp+1:time_arrival_samp+arrival_dur_ms/1000*fs,:) = arrival_window;

srir_arrival_removed = srir_in .* (1-srir_extract_window);
% question of whether to use srir_in or srir_unAlt. Unalt can produce
% build up of reflections
srir_arrival = srir_in(time_arrival_samp+1:time_arrival_samp+arrival_dur_ms/1000*fs,:) .* arrival_window;
% srir_arrival = srir_unAlt(time_arrival_samp+1:time_arrival_samp+arrival_dur_ms/1000*fs,:) .* arrival_window;

end


function[doa_est,normalized_doa_est_P_dB,P_pwd] = plot_doa_horiz(srir,fs)
% plot horizontal doa

% Tuneable parameters:
res_deg_azi = 2;
res_deg_ele = 90;
order = sqrt(size(srir,2))-1;
nSrc = 7;
numSamps = fs*0.3; % 0.3 seconds
highPassFilterFreq = 2000;
kappa = 40;

grid_dirs = grid2dirs(res_deg_azi,res_deg_ele,0,0); % Grid of directions to evaluate DoA estimation

% remove elevated directions -- only for a res_deg_ele of 90
grid_dirs = grid_dirs(2:end-1,:);

P_src = diag(ones(numSamps,1));
[~,filtHi,~] = ambisonic_crossover(highPassFilterFreq,fs);

doa_est = zeros(nSrc,2,length(srir(1,1,:)));
doa_est_P = zeros(nSrc,length(srir(1,1,:)));
P_pwd = zeros(size(grid_dirs,1),length(srir(1,1,:)));

for i = 1:size(srir,3)
    Y_src = filter(filtHi,1,srir(1:numSamps,:,i)); % high pass filter
    stVec = Y_src';
    sphCOV = stVec*P_src*stVec' + 1*eye((order+1)^2)/(4*pi);

    % DoA estimation
    [P_pwd(:,i), est_dirs_pwd,est_dirs_P] = sphPWDmap(sphCOV, grid_dirs, nSrc,kappa);
    est_dirs_pwd = est_dirs_pwd*180/pi; % convert to degs from rads

    % flip -ve values near -180 to +ve values
    negativeFlipLimit = -178;
    for j = 1:length(est_dirs_pwd(:,1))
        for k = 1:length(est_dirs_pwd(1,:))
            if est_dirs_pwd(j,k) < negativeFlipLimit
                est_dirs_pwd(j,k) = est_dirs_pwd(j,k) + 360;
            end
        end
    end
    doa_est(:,:,i) = est_dirs_pwd;
    doa_est_P(:,i) = est_dirs_P;
end

% normalized_doa_est_P = doa_est_P ./ max(doa_est_P,[],1); % NORMALISE
normalized_doa_est_P = doa_est_P;% DONT NORMALISE
normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;

% plot_thresh = -20; % for normalised powers
plot_thresh = min(normalized_doa_est_P_dB(:)); % for unnormalised powers

normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;
cmap = flip(parula(256));
c_truncation = 20; c = cmap(c_truncation:end,:); % truncate yellow

doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);
cspace = linspace(0,1,size(c,1));
for i = 1:size(srir,3)
    doa_color(:,:,i) = interp1(cspace, c(:,i), doa_01);
end

hold off
for i = nSrc:-1:1
    s = scatter(squeeze(doa_est(i,1,:)),1:size(srir,3),150,...
        squeeze(doa_color(i,:,:)),"x",'LineWidth',3);
    hold on
end

xlabel('Azimuth (°)');ylabel('');
xlim([negativeFlipLimit (negativeFlipLimit+360)]);xticks(-180:45:180);
colormap(c);
k = colorbar;
xlabel(k,'Normalised power (dB)');caxis([plot_thresh 0]);
set(gcf, 'Color', 'w');pbaspect([1.7 1 1]);
box on;grid on;set(gca,'FontSize',16)

ylim([0.5 size(srir,3)+0.5]);
yticks(1:size(srir,3));
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
end

function[doa_deg] = doa_iv(srir) % gives answer in degrees
% intensity vector DoA analysis. (x channel times [x, y, z] channel). Gives
% the amount of energy in each axis. If it's pointing forwards, you get
% positive energy in the X axis. if it's pointing left, you get positive
% energy in the y axis. up -> positive energy in z axis
doa_samp = srir(:, 1) .* [srir(:, 4),...
    srir(:, 2), srir(:, 3)];
idxNonZero = find(sum(doa_samp ~= 0, 2));
if isempty(idxNonZero)
    doa = mean(doa_samp(:, :),1);
else
    doa = mean(doa_samp(idxNonZero, :),1);
    doa = doa ./ vecnorm(doa, 2, 2);
end
[doa1(1,1),doa1(2,1),~] = cart2sph(doa(1),doa(2),doa(3));
doa_deg = rad2deg(doa1);
end

function [P_pwd,doa_est,doa_est_P,grid_dirs] = get_pwd(srir,fs)

degreeResolution = 2;
order = sqrt(size(srir,2))-1;
nSrc = 7;
numSamps = fs*0.3; % 0.3 seconds
highPassFilterFreq = 3000;
kappa = 40;

grid_dirs = grid2dirs(degreeResolution,degreeResolution,0,0); % Grid of directions to evaluate DoA estimation
P_src = diag(ones(numSamps,1));
[~,filtHi,~] = ambisonic_crossover(highPassFilterFreq,fs);

doa_est = zeros(nSrc,2);
doa_est_P = zeros(nSrc,1);
P_pwd = zeros(size(grid_dirs,1),1);

Y_src = filter(filtHi,1,srir(1:numSamps,:)); % high pass filter
stVec = Y_src';
sphCOV = stVec*P_src*stVec' + 1*eye((order+1)^2)/(4*pi);

% DoA estimation
[P_pwd(:), est_dirs_pwd,est_dirs_P] = sphPWDmap(sphCOV, grid_dirs, nSrc,kappa);

est_dirs_pwd = est_dirs_pwd*180/pi; % convert to degs from rads
% flip -ve values near -180 to +ve values
negativeFlipLimit = -170;
for j = 1:length(est_dirs_pwd(:,1))
    for k = 1:length(est_dirs_pwd(1,:))
        if est_dirs_pwd(j,k) < negativeFlipLimit
            est_dirs_pwd(j,k) = est_dirs_pwd(j,k) + 360;
        end
    end
end
doa_est(:,:) = est_dirs_pwd;
doa_est_P(:) = est_dirs_P;
end