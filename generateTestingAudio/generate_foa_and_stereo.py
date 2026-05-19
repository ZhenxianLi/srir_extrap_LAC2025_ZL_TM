# -*- coding: utf-8 -*-
"""
此脚本用于生成一个一级 Ambisonics (FOA) B-Format 音频文件，
并将其解码为一个立体声音频文件以供验证。

此版本已修改为使用内置的 `wave` 模块，不再需要 `scipy` 库。
"""

import numpy as np
import wave
import os


# --- 1. 参数设置 ---
SAMPLE_RATE = 48000  # 采样率 (Hz)
DURATION = 5         # 持续时间 (秒)
FREQUENCY = 440      # 正弦波频率 (Hz, A4音)
AZIMUTH_DEG = 90     # 方位角 (度), 0°为正前方, 90°为正左方, -90°为正右方
ELEVATION_DEG = 0    # 俯仰角 (度), 0°为水平面

OUTPUT_FILENAME_FOA = "FOA_Sine_90L_0U.wav"
OUTPUT_FILENAME_STEREO = "FOA_Decoded_Stereo.wav"
AMPLITUDE = 0.7      # 幅度，留出一些动态余量


def write_wav(filename, rate, data_int16):
    """使用内置 wave 模块将 numpy int16 数组写入 WAV 文件。"""
    with wave.open(filename, 'w') as wf:
        if len(data_int16.shape) == 1:
            num_channels = 1
        else:
            num_channels = data_int16.shape[1]

        wf.setnchannels(num_channels)
        wf.setsampwidth(2)  # 16-bit -> 2 bytes
        wf.setframerate(rate)
        wf.writeframes(data_int16.tobytes())

print("开始生成音频信号...")

# --- 2. 生成单声道声源信号 ---
num_samples = int(SAMPLE_RATE * DURATION)
t = np.linspace(0., DURATION, num_samples, endpoint=False)
source_signal = AMPLITUDE * np.sin(2 * np.pi * FREQUENCY * t)
print(f"已生成 {DURATION} 秒, {FREQUENCY} Hz 的正弦波声源。")

# --- 3. 将声源编码为一级 Ambisonics (FOA) ---
azimuth_rad = np.deg2rad(AZIMUTH_DEG)
elevation_rad = np.deg2rad(ELEVATION_DEG)
w_channel = source_signal
x_channel = source_signal * np.cos(azimuth_rad) * np.cos(elevation_rad)
y_channel = source_signal * np.sin(azimuth_rad) * np.cos(elevation_rad)
z_channel = source_signal * np.sin(elevation_rad)
foa_signal = np.stack([w_channel, x_channel, y_channel, z_channel], axis=-1)
print(f"已将声源编码为FOA格式，声源位置: 方位角={AZIMUTH_DEG}°, 俯仰角={ELEVATION_DEG}°")

# --- 4. 保存为 4 通道 .wav 文件 ---
foa_signal_int16 = np.int16(foa_signal * 32767)
write_wav(OUTPUT_FILENAME_FOA, SAMPLE_RATE, foa_signal_int16)
print(f"成功保存4通道FOA文件: {os.path.abspath(OUTPUT_FILENAME_FOA)}")

# --- 5. 验证：将 FOA 解码为立体声 ---
print("\n开始解码以进行验证...")
left_channel = 0.5 * (foa_signal[:, 0] + foa_signal[:, 2])
right_channel = 0.5 * (foa_signal[:, 0] - foa_signal[:, 2])
stereo_signal = np.stack([left_channel, right_channel], axis=-1)

# --- 6. 保存为立体声 .wav 文件 ---
stereo_signal_int16 = np.int16(stereo_signal * 32767)
write_wav(OUTPUT_FILENAME_STEREO, SAMPLE_RATE, stereo_signal_int16)
print(f"成功保存解码后的立体声文件: {os.path.abspath(OUTPUT_FILENAME_STEREO)}")
print("\n验证完成！")
print(f"请播放 '{OUTPUT_FILENAME_STEREO}' 文件，您应该会听到声音完全来自左声道。")