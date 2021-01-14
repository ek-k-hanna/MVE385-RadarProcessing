function data = generate_clutter(amp, freq, angle, angle_bw, ...
                                    samples, antennas, sampling_freq)    
    % Sampling offset for antenna phase + 
    % frequency modulation of the antenna angle.
    wave_length = sampling_freq/freq;
    angle_indexes_offset = ((0:antennas - 1)*wave_length/2*angle/90) + ...
                            sin(pi/2 + pi*(0:antennas - 1)/samples)*wave_length/2*angle_bw/2/90;
    
    % Sample points with considerations of antenna phase + clutter BW
    x = repmat((0:samples - 1)', 1, antennas);
    x = bsxfun(@plus, angle_indexes_offset, x);
    
    amp = 10^(amp/20);
    
    % Wave equation
    data = amp*exp(1i*2*pi*x*freq/sampling_freq);
end
