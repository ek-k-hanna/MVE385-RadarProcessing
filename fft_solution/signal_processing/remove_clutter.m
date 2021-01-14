function S = remove_clutter(S, cut_off_freq, PRF, attenuation, slope)
    for i = 1:size(S,2)
        S(:,i) = highpass(S(:,i), cut_off_freq, PRF, 'StopbandAttenuation', attenuation, 'Steepness', slope);
    end
end