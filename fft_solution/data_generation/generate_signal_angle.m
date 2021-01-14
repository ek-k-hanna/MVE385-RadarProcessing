function data = generate_signal_angle(amp, freq, angle, samples, antennas, sampling_freq)
    % Sampling offset for antenna phase
    wave_length = sampling_freq/freq;
    angle_indexes_offset = ((0:antennas - 1)*wave_length/2*angle/90);
    
    % Sample points with considerations of antenna phase
    x = repmat((0:samples - 1)', 1, antennas);
    x = bsxfun(@plus, angle_indexes_offset, x);
    
    amp = 10^(amp/20);
    
    % Wave equation
    data = amp*exp(1i*2*pi*x*freq/sampling_freq);
end