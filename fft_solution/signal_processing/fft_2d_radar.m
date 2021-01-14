function fft_data_2d = fft_2d_radar(S)
    fft_data_1d = fftshift(fft(S, [], 1), 1);
    fft_data_2d = fftshift(fft(fft_data_1d, [], 2), 2);
    %fft_data_2d = abs(fftshift(fft(fft(S')'), 1));
end