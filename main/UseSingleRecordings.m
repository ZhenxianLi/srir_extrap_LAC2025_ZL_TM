% adjust from RotateIR.m
% scale all Early reflect
% from single recording
%--------------------------------------------------------------------------------------------------------------------------
%read and plot sofa file

clc;
clear all;
close all;

% Start SOFA and load function library

% get the now path
parentFolder = fileparts(pwd);

addpath(fullfile(parentFolder, 'API_MO'));
addpath(fullfile(parentFolder, 'Spherical-Harmonic-Transform-master'));
addpath(fullfile(parentFolder, 'Higher-Order-Ambisonics-master'));
addpath(genpath(fullfile(parentFolder, 'binaural_ambisonic_preprocessing-main')));
addpath(fullfile(parentFolder, 'audio_samples'));
SOFAstart;


%run this to load decoder
if ~exist('aio_flag')
    load_ambisonic_configuration
end

% Load Variable Acoustics 6DoF Dataset (most reverberant set)

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
%SOFAplotGeometry(sofa1);   % plot the source and listen position

%% check inf
sofa1.GLOBAL_RoomType;
sofa1.GLOBAL_SOFAConventionsVersion;
sofa1.ListenerPosition;
sofa1.SourcePosition;
sofa1.EmitterPosition;
size(sofa1.Data);

%cheak the first  by hand
firstReflectionIndex=[487;697;674;528;579;621;606;515;616;726;499;658;861;543;742;939;971;1146;1071;900;1171];
%% calculate the scale parameter & sound speed
sizeIR=size(sofa1.Data.IR);
N=sizeIR(1);    
N=15;     %only facing=15

IROrder=sizeIR(2);
L=sizeIR(3);

% find direct sound and fisrt early ref sound arrive time

speedSum=0;
distance=zeros(N,1);
directSoundValure=zeros(N,1);
for n=1:N
    IR_n=squeeze(sofa1.Data.IR(n,1,:));% n,1st order,full length
    [directSoundIndex,directSoundValure(n),pks,lcs]=findDirectSound(IR_n);
    distance(n)=sqrt(sum((sofa1.ListenerPosition(n,:)-sofa1.SourcePosition(n,:)).^2));
    directSoundTime(n)=directSoundIndex;     %time of direct sound(in samples)
    speedSum=(speedSum+distance(n)/directSoundTime(n));
    speed=speedSum/n;% in samples

end

%plot all N sound
pic1=0;
if pic1
    for n=1:N
        figure;
        subplot(7,3,n);%3 speaker; 7 mic
        plot(abs(IR_n));
        hold on;
        xline(directSoundIndex,'--g');
        xline(firstReflectionIndex(n),'--g');
        ylim([0,0.5]);
        xlim([0,4000]);
        hold on;
        plot(lcs,pks);
        title(num2str(n))

    end


    %end for loop
end

% generate new IR-------------------------------
%% cut and scale
% start from 9, S3L3 ,center of room
% which postion use to

% generate*************************************************************************************
%target position is the Speaker:usespeaker ;Mic:uselistener
useSpeaker=1;%1-3
useListener=5;%1-7

%*********************************************************************************************

controlRecordNum=useListener*3-3+useSpeaker;%target
directSoundCutLeft=50;
directSoundCutRight=300;
earlyRefCutLength=400;
%read IR and postion
IR_original=sofa1.Data.IR(9,:,:);

SourcerPoint_Origin=sofa1.SourcePosition(9,:);
ListenerPoint_Origin=sofa1.ListenerPosition(9,:);

SourcerPoint_generate=sofa1.SourcePosition(controlRecordNum,:);
ListenerPoint_generate=sofa1.ListenerPosition(controlRecordNum,:);% the controlRecordNum recording it target
arrival_time_original=directSoundTime(9);
%Am_original=directSoundValure(9);
% cut
directSound=IR_original(:,:,arrival_time_original-directSoundCutLeft:arrival_time_original+directSoundCutRight);

%fit distance - amiptude
%use rms
Nface=15;%only facing 15
rmsDirectSoundValure=zeros(Nface,1);

for n=1:Nface
    readIR_forFit=sofa1.Data.IR(n,:,:);
    arrival_time_forFit=firstReflectionIndex(n);
    ds=readIR_forFit(:,1,arrival_time_forFit-directSoundCutLeft:arrival_time_forFit+directSoundCutRight);
    rmsDirectSoundValure(n)=rms(squeeze(ds));
end
fun = @(x)sseval(x,distance(1:Nface,1),rmsDirectSoundValure);
x0 = rand(2,1);
bestx = fminsearch(fun,x0);

pic2=1;
if pic2
    figure;% check the fit
    scatter(distance,rmsDirectSoundValure);
    for n=1:N
        text(distance(n),rmsDirectSoundValure(n),num2str(n));
    end
    hold on;
    dislist=0:0.2:10;
    plot(dislist,bestx(1)*exp(-(bestx(2))*dislist));
    xlabel('distace/m')
    ylabel('rms Amplitude of Direct Sound')

end

%compute the time and Amiplitude
distance_generate=sqrt(sum((SourcerPoint_generate-ListenerPoint_generate).^2));
distance_origianl=sqrt(sum((SourcerPoint_Origin-ListenerPoint_Origin).^2));
arrival_time_generate=floor(distance_generate/speed);       %  in samples
Am_generate=bestx(1)*exp(-(bestx(2))*distance_generate);
Am_original=bestx(1)*exp(-(bestx(2))*distance_origianl);
directSound=directSound*(Am_generate/Am_original); %scale

% re-join
A=zeros(25,arrival_time_generate-directSoundCutLeft-1);
B=squeeze(directSound);         %direct sound
C=squeeze(IR_original(:,:,arrival_time_generate+directSoundCutRight+1:arrival_time_generate+directSoundCutRight+earlyRefCutLength));%early ref
%C=C*(Am_generate/Am_original)/0.25;%scale C
D=squeeze(IR_original(:,:,arrival_time_generate+directSoundCutRight+earlyRefCutLength+1:end));%reverb

IR_scale=[A,B,C,D];

%% plot and listen
%plot after scale
pic3=0;
if pic3
    %sound(IR_scale(1,:),fs);
    figure;

    subplot(3,1,1);
    plot(abs(squeeze(IR_original(:,1,:))));
    legend('IR original');
    xlim([0 3000]);
    ylim([0 0.5]);

    subplot(3,1,2);
    plot(abs(squeeze(IR_scale(1,:))));
    legend('IR generate');
    xlim([0 3000]);
    ylim([0 0.5]);

    subplot(3,1,3);
    plot(abs(squeeze(sofa1.Data.IR(controlRecordNum,1,:))));
    legend('IR record');
    xlim([0 3000]);
    ylim([0 0.5]);
end
%% rotate

% compute the sph
xyz1=SourcerPoint_Origin - ListenerPoint_Origin;
xyz2=SourcerPoint_generate-ListenerPoint_generate;
[az1,el1,r1]=cart2sph(xyz1(1),xyz1(2),xyz1(3));
[az2,el2,r2]=cart2sph(xyz2(1),xyz2(2),xyz2(3));
yaw=(az2-az1)/pi*180;
pitch=(el2-el1)/pi*180;
roll=0;

%reshape hoasig
hoasig=B.'; %B- reshape direct sound

% rotate
hoasig_rot = rotateHOA_N3D(hoasig, yaw, pitch, roll);
B_rot=hoasig_rot.';

%re-join
IR_rot=[A,B_rot,C,D];


%% plot and listen
%plot after rotate
pic3=1;
if pic3
    %sound(IR_scale(1,:),fs);
    figure;

    subplot(3,1,1);
    plot(abs(squeeze(IR_original(:,1,:))));
    legend(['IR original ',num2str(3),num2str(3)]);
    xlim([0 8000]);
    ylim([0 0.5]);

    subplot(3,1,2);
    plot(abs(squeeze(IR_rot(1,:))));
    legend('IR rot');
    xlim([0 8000]);
    ylim([0 0.5]);

    subplot(3,1,3);
    plot(abs(squeeze(sofa1.Data.IR(controlRecordNum,1,:))));
    legend(['IR record ',num2str(useSpeaker),num2str(useListener)]);
    xlim([0 8000]);
    ylim([0 0.5]);
end
%% render binaural and play
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
    IR_Record=squeeze(sofa1.Data.IR(controlRecordNum,:,:)) ;
    IR_original=squeeze(sofa1.Data.IR(9,:,:));

    binIR_original=binSound(IR_original,SH_ambisonic_binaural_decoder);
    binIR_generate=binSound(IR_rot,SH_ambisonic_binaural_decoder);
    binIR_record=binSound(IR_Record,SH_ambisonic_binaural_decoder);
    gaptime=zeros(0.3*fs,2);
    soundsc([binIR_original;gaptime;binIR_generate;gaptime;binIR_record],Fs);
    soundsc([binIR_generate;gaptime;binIR_record],Fs);% only last 2 sound

end

%play controlRecord then generate IR conv with dry guitar
% 1.S3L3 2.generete 3.target recoeding
doBinRenderSONG=0;
if doBinRenderSONG
    %song_Dry=audioread('speechdirectsound_48.wav');
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


    audiowrite(['S',num2str(useSpeaker),'L',num2str(useListener),'_Original33','.wav'],binSong_originalIR/Max,fs);
    audiowrite(['S',num2str(useSpeaker),'L',num2str(useListener),'_Rotate','.wav'],binSong_generateIR/Max,fs);
    %audiowrite(['S',num2str(useSpeaker),'L',num2str(useListener),'_Record','.wav'],binSong_recordIR/Max,fs);

end


%% declare function
%Find  direct sound loc and val
function [locD,ValD,pks,lcs]=findDirectSound(ir)
absir=abs(ir);
noiseValuse=max(absir);
[pks,lcs]=findpeaks(absir,"MinPeakDistance",10,MinPeakHeight=0.2*noiseValuse);  %denosing,find peak
% Find the index of the peak representing the direct sound
locD=lcs(1);
ValD=pks(1);
end

% find A0 and alpha for distance - amipitude fit

function sse = sseval(x,tdata,ydata)
A = x(1);
lambda = x(2);
sse = sum((ydata - A*exp(-lambda*tdata)).^2);
end


function binIR=binSound(IR_test,SH_ambisonic_binaural_decoder)% render to a bin audio and play it
ambisonic_soundscape=IR_test.';

binaural_ambisonic_render = zeros(length(SH_ambisonic_binaural_decoder(1,:,1))+length(ambisonic_soundscape)-1,length(SH_ambisonic_binaural_decoder(1,1,:)));

% convolve each channel of the encoded signal with the decoder signal and sum the result
for i = 1:length(SH_ambisonic_binaural_decoder(:,1,1))

    binaural_ambisonic_render(:,1) = binaural_ambisonic_render(:,1) +  conv(SH_ambisonic_binaural_decoder(i,:,1),ambisonic_soundscape(:,i) );
    binaural_ambisonic_render(:,2) = binaural_ambisonic_render(:,2) +  conv(SH_ambisonic_binaural_decoder(i,:,2),ambisonic_soundscape(:,i) );
end

% compare the two binaural decoders by listening to both binaural renders consecutively
binIR=binaural_ambisonic_render;
end





