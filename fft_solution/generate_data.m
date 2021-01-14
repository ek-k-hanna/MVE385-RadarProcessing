clear;

addpath(genpath('data_generation'))

% DATA CONSTANT SETTINGS
samples = 64;
antennas = 24;
n = samples*antennas;
sampling_freq = 1000;

% DATA VARIABLES

% TARGET 1
target1_amp_dB = 0;
target1_freq = 250; % +- 500Hz
target1_angle = 45; % +- 90*

% TARGET 2
multiple_targets_flag = true;
target2_amp_dB = 0;
target2_freq = -125; % +- 500Hz
target2_angle = -22.5; % +- 90*

% NOISE
noise_amp_per_bin = true; % Else energy for whole BW
noise_amp_full_BW_dB = 23;
noise_amp_bin_dB = -20;

% CLUTTER
clutter_flag = false;
clutter_amp_full_BW_dB = 20; % Energy get smudged over angle BW
clutter_freq = 1000/64; % 1000/64 to 1000/64*3 Hz is ok.
clutter_angle_centre = 90/12*6; % +- 90
clutter_angle_bw = (180/24)*180; % BLACK MAGIC - MAX ??? (180/24)*180

% FLAG MODIFIERS
if ~multiple_targets_flag
    target2_amp_dB = -inf;
end

if noise_amp_per_bin
    noise_amp_full_BW_dB = noise_amp_bin_dB + 20*log10(n);
else
    noise_amp_bin_dB = noise_amp_full_BW_dB - 20*log10(n);
end

if ~clutter_flag
    clutter_amp_full_BW_dB = -inf;
end

% GENERATE DATA
target1 = generate_signal_angle(target1_amp_dB, target1_freq, target1_angle, samples, antennas, sampling_freq);
target2 = generate_signal_angle(target2_amp_dB, target2_freq, target2_angle, samples, antennas, sampling_freq);
noise = generate_noise(noise_amp_full_BW_dB, samples, antennas);
clutter = generate_clutter(clutter_amp_full_BW_dB, clutter_freq, clutter_angle_centre, clutter_angle_bw, ...
                            samples, antennas, sampling_freq);
data = target1 + target2 + noise + clutter;

% Check amplitudes in time domain
target1_power = time_domain_estimate_power(target1)
target2_power = time_domain_estimate_power(target2)
noise_power = time_domain_estimate_power(noise)
noise_power_per_bin = noise_power - 20*log10(n)
clutter_power = time_domain_estimate_power(clutter)
SNR_greatest_target = max(target1_power - noise_power_per_bin, target2_power - noise_power_per_bin)

% Check amplitudes in 1D freq domain
% WARNING, dont get 1D to work.
%S1 = fftshift(fft(target1, [], 1), 1);
%10*log10(sum(abs(S1(:,1)).^2)/64^2)

%S1 = fftshift(fft(noise, [], 1), 1);
%20*log10(sum(abs(S1(:,1)).^2)/64^2/n)

% Check amplitudes in 2D freq domain
S1 = fftshift(fft(target1, [], 1), 1);
S2 = fftshift(fft(S1, [], 2), 2);
20*log10(sum(abs(S2).^2, 'all')/64^2/24^2)

S1 = fftshift(fft(noise, [], 1), 1);
S2 = fftshift(fft(S1, [], 2), 2);
20*log10(mean(abs(S2).^2, 'all')/64.^2/24.^2)

% FFT
S1 = fftshift(fft(data, [], 1), 1);
S2 = fftshift(fft(S1, [], 2), 2);
S = abs(S2);



%% 1D FFT DISPLAY
plot(20*log10(abs(S1(:,1))))

%% 2D FFT DISPLAY
imagesc(S);
xlabel('Angle [bins]')
ylabel('Doppler towards the antenna [bins]')
colorbar;

%% Compute SNR (Only works for single target with no clutter)

if (~multiple_targets_flag & ~clutter_flag)
    [M, I_target] = max(S, [], 'all', 'linear');
    I_noise = setdiff(1:numel(S), I_target);

    target_amp = abs(S(I_target)).^2
    noise_total_amp = sum(abs(S(I_noise)).^2, 'all')/(n-1)*n
    noise_bin_amp = mean(abs(S(I_noise)).^2, 'all')

    SNR_total = 20*log10(target_amp/noise_total_amp)
    SNR_bin = 20*log10(target_amp/noise_bin_amp)
else
    disp('WARNING: Cant compute SNR for multiple targets or clutter.')    
end

%% FUNCTIONS

function power = time_domain_estimate_power(data)
    power = mean(abs(data).^2, 'all');
    power = 20*log10(power);
end







