function S_windowed = window_2d(S, window_type)
    n_samples_row = size(S, 1);
    w_row = window_type(n_samples_row);
    window_offset_row = n_samples_row/sum(w_row);
    
    n_samples_col = size(S, 2);
    w_col = window_type(n_samples_col);
    window_offset_col = n_samples_col/sum(w_col);
    
    S_windowed = S;
    for i = 1:size(S,1) %% doppler | time axis
        for j = 1:size(S,2)
            S_windowed(i,j) = S(i,j).*w_row(i,1)*w_col(j,1);
        end
    end
    S_windowed = S_windowed*window_offset_row*window_offset_col;
end
