function display_fft_bins_2d(fft_data)
    imagesc(fft_data);
    xlabel('Angle [bins]')
    ylabel('Doppler towards the antenna [bins]')
    colorbar;
end
