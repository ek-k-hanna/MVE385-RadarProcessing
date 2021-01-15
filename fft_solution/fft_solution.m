clear

%rng('default')
addpath(genpath('data_generation'))
addpath(genpath('signal_processing'))
addpath(genpath('improvement_methods'))
addpath(genpath('visualization'))

% SYSTEM PARAMETERS
carrier_freq = 1.3E9;
PRF = 1E3;
speed_of_light = 3E8;

% DATA SETTINGS
generate_data = true;

% FFT ALGORITHMS SETTINGS
window_flag = true;

% VISUALIZATION
visualization_flag = true;

% DATA CONSTANT SETTINGS
samples = 64;
antennas = 24;
n = samples*antennas;
sampling_freq = 1000;

%% LOAD DATA

% SAAB data
if ~generate_data

% Task 1
%raw_data = load('data/Example_1a.mat');
%raw_data = load('data/Example_1b.mat');
%raw_data = load('data/Example_1c.mat');
%raw_data = load('data/Example_1d.mat');
% Task 2
raw_data = load('data/Example_2a_target_in_clutter.mat');
%raw_data = load('data/Example_2b_target_in_clutter.mat');
% Task 3
%raw_data = load('data/Example_3_two_targets.mat');
S = raw_data.S;
clear raw_data;

else
% GEN DATA

% TARGETS
linear_freq = true; % Otherwise velocity and az angle.

% TARGET 1
target1_amp_dB = 0;
if linear_freq
    target1_freq = 500/32*6.3; % +-500 != 0.
    target1_vel = doppler_freq_to_velocity(target1_freq, carrier_freq);
    target1_angle = 90/12*5.2; % +-90
else
    velocity = 15; % +-57.6923 != 0.
    target1_freq = velocity_to_doppler_freq(velocity, carrier_freq);
    az_deg = 45; % +-90
    target1_angle = az_deg_to_angle_freq(az_deg);
end

% TARGET 2
multiple_targets_flag = false;
target2_amp_dB = 0;
if linear_freq
    target2_freq = 500/32*5; % +-500 != 0.
    target2_angle = 90/12*7; % +-90
else
    velocity = -30; % +-57.6923 != 0.
    target2_freq = velocity_to_doppler_freq(velocity, carrier_freq);
    az_deg = -70; % +-90
    target2_angle = az_deg_to_angle_freq(az_deg);
end

if (target1_freq == 0) || (target2_freq == 0)
   disp('WARNING: Velocity can not be zero!')   
end

% NOISE
noise_amp_per_bin = true; % Else energy for whole BW
noise_amp_full_BW_dB = 23;
noise_amp_bin_dB = -50;

% CLUTTER
clutter_flag = false;
clutter_amp_full_BW_dB = 10; % Energy get smudged over angle BW
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
    clutter_amp_full_BW_dB = -10000;
end

% GENERATE DATA
target1 = generate_signal_angle(target1_amp_dB, target1_freq, target1_angle, samples, antennas, sampling_freq);
target2 = generate_signal_angle(target2_amp_dB, target2_freq, target2_angle, samples, antennas, sampling_freq);
noise = generate_noise(noise_amp_full_BW_dB, samples, antennas);
clutter = generate_clutter(clutter_amp_full_BW_dB, clutter_freq, clutter_angle_centre, clutter_angle_bw, ...
                            samples, antennas, sampling_freq);
S = target1 + target2 + noise + clutter;
end

% DISPLAY ORIGINAL DATA

abs_fft_data = abs(fft_2d_radar(S));
if visualization_flag
    figure;
    display_fft_bins_2d(abs_fft_data);
end

%% Signal Processing

% FILTER OUT CLUTTER
cut_off_freq = PRF/128*3;
attenuation = 45;
slope = 0.5;
S_hp = remove_clutter(S, cut_off_freq, PRF, attenuation, slope);

% WINDOW
if window_flag
    S_windowed2 = window_2d(S_hp, @hamming);
else
    S_windowed2 = S_hp;
end
    
% FFT
fft_data = fft_2d_radar(S_windowed2);

% DISPLAY FILTERED DATA
if visualization_flag
    figure;
    abs_fft_data = abs(fft_data);
    display_fft_bins_2d(abs_fft_data);
end

%% Kernel

edge_detection_kernel_3x3 = [0 -1 0; -1 4 -1; 0 -1 0];
%edge_detection_kernel_3x3 = [-1 -1 -1; -1 8 -1; -1 -1 -1];

C = conv2(abs_fft_data, edge_detection_kernel_3x3, 'same');

indexes = 1:numel(C);
clutter_indexes = (32:34)' + size(C,1).*((1:size(C,2))-1);
none_clutter_indexes = setdiff(indexes, clutter_indexes);

% STATISTICS
C_mean = mean(C(none_clutter_indexes));
C_std = std(C(none_clutter_indexes));

C_std_multipler = 8; % 4 on SAABs data. 8 on Chalmers groups synth.
signal_indexes = find(C > C_mean + C_std * C_std_multipler); % MAGIC NUMBER 4 WORKS?
signal_indexes = setdiff(signal_indexes, clutter_indexes);
signal_indexes_expanded = [signal_indexes; ...
                            signal_indexes - 1; signal_indexes + 1; ...
                            signal_indexes - 64; signal_indexes + samples];
signal_indexes_expanded = intersect(signal_indexes_expanded, 1:(samples*antennas));
noise_indexes = setdiff(indexes, union(signal_indexes_expanded, clutter_indexes));

% COMPUTE NOISE LEVEL
noise = mean(abs_fft_data(noise_indexes));
noise_power = mean(abs(fft_data(noise_indexes)).^2, 'all');

% DISPLAY KERNEL IMAGE
if visualization_flag
    figure;
    display_fft_bins_2d(C);
end

% DISPLAY NOISE
if visualization_flag
    C = abs_fft_data;
    C(noise_indexes) = 1;
    C(setdiff(indexes, noise_indexes)) = 0;
    figure;
    display_fft_bins_2d(C);
end

% DISPLAY SIGNAL
if visualization_flag
    C = abs_fft_data;
    C(signal_indexes) = 1;
    C(setdiff(indexes, signal_indexes)) = 0;
    figure;
    display_fft_bins_2d(C);
end

%% Find peaks

fft_data_padded = pad_matrix(fft_data, noise);
C = abs(fft_data_padded);
signal_indexes_paded = pad_index(signal_indexes, samples);

if isempty(signal_indexes_paded)
    disp('WARNING: Did not find any target!');
else

    signals = zeros(size(signal_indexes))';
    j = 1;
    for i = signal_indexes_paded'
        if all(C(i) > C(i + [-1; 1; -(samples + 2); (samples + 2)]))
            signals(j) = 1;
        end
        j = j + 1;
    end

    targets = signal_indexes(logical(signals));

    % WHERE ARE THE PEAKS?
    if visualization_flag
        C_peaks = abs_fft_data;
        C_peaks(targets) = 1;
        C_peaks(setdiff(indexes, targets)) = 0;
        figure;
        display_fft_bins_2d(C_peaks);
    end

    %% Improving
    for i = 1:numel(targets)
        disp(['Computing parameters for target index: ', i + '0']); % ASCII magic

        %find coordinates of detected item in matrix
        [rIdx, cIdx] = get_row_col(pad_index(targets(i), samples), samples + 2);

        lagrange_1d_pointMax = lagrange_1d(C, rIdx, cIdx);
        lagrange_2d_pointMax = lagrange_2d(C, rIdx, cIdx);
        weighted_average_pointMax = weighted_average(C, rIdx, cIdx);

        signal_amp = abs(C(rIdx, cIdx)) - noise;

        signal_log20db = 20*log10(signal_amp.^2/samples^2/antennas^2);
        noise_log20db = 20*log10(noise_power/samples^2/antennas^2);
        meassured_SNR = signal_log20db - noise_log20db;

        rIdx = rIdx - 1;
        cIdx = cIdx - 1;
        lagrange_1d_pointMax = lagrange_1d_pointMax - 1;
        lagrange_2d_pointMax = lagrange_2d_pointMax - 1;
        weighted_average_pointMax = weighted_average_pointMax - 1;

        % Convert to freq, then convert to velocity and angle
        vel_freq_centre_bin = bin_to_freq(rIdx, samples, PRF)
        ang_freq_centre_bin = bin_to_freq(cIdx, antennas, 180)
        
        vel_freq_lagrange_1d = bin_to_freq(lagrange_1d_pointMax(1), samples, PRF)
        ang_freq_lagrange_1d = bin_to_freq(lagrange_1d_pointMax(2), antennas, 180)

        vel_freq_lagrange_2d = bin_to_freq(lagrange_2d_pointMax(1), samples, PRF)
        ang_freq_lagrange_2d = bin_to_freq(lagrange_2d_pointMax(2), antennas, 180)

        vel_freq_weighted_average = bin_to_freq(weighted_average_pointMax(1), samples, PRF)
        ang_freq_weighted_average = bin_to_freq(weighted_average_pointMax(2), antennas, 180)

        %meassured_vel = doppler_freq_to_velocity(vel_freq, carrier_freq)
        %meassured_ang = angle_freq_to_az_deg(ang_freq)
    end
end

%% Close open windows
close all

%% FUNCTIONS

% DATA REPRESENTATION CONVERSION

function [row, col] = get_row_col(index, rows)
    row = mod(index, rows);
    col = (index - row)/rows + 1;
end

% PADDING WITH NOISE INDEXING CONVERSION

function index = pad_index(index, rows)
    mod_index = mod(index - 1, rows) + 1;
    index = rows + 2 + ((index - mod_index)/rows)*(rows+2) + 1 + mod_index;
end

function index = unpad_index(index, rows)
    mod_index = mod(index, rows);
    index = index - 1 - (((index - mod_index)/rows)-1)*2 - rows;
end

function M = pad_matrix(S, value)
    M = ones(size(S) + 2) * value;
    M(2:(size(S,1)+1), 2:(size(S,2)+1)) = S;
end

function M = unpad_matrix(S)
    M = S(2:(size(S,1)-1), 2:(size(S,2)-1));
end

% BIN AND FREQ CONVERSION

function freq = bin_to_freq(bin, bin_bw, freq_bw)
    freq = (bin - (bin_bw/2 + 1))*freq_bw/bin_bw;
end

% FREQ AND PARAMETER CONVERSION

function angle_freq = az_deg_to_angle_freq(az_deg)
    angle_freq = sin((az_deg)/180*pi)*90;
end

function az_deg = angle_freq_to_az_deg(angle_freq)
    az_deg = asin(angle_freq/90)*180/pi;
end

function doppler_freq = velocity_to_doppler_freq(velocity, carrier_freq)
    speed_of_light = 3E8;
    doppler_freq = carrier_freq * velocity / speed_of_light * 2;
end

function velocity = doppler_freq_to_velocity(doppler_freq, carrier_freq)
    speed_of_light = 3E8;
    velocity = doppler_freq * speed_of_light / carrier_freq / 2;
end

