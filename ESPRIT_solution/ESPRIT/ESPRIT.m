function frequencies = ESPRIT(signal,number_frequecies,signal_length)
    
    L = floor((signal_length+1)/2); % optimal separation condition 
    K = signal_length - L + 1; 
    % create hankell matrix and sub-matrices 
    H = hankel(signal);
    H = H(1:L+1,1:K);
    H1 = H(1:L,:);
    H2 = H(2:L+1,:);
    
    % SVD 
    [U,Sigma,V] = svd(H1);
    U1 = U(:,1:number_frequecies);  % truncate
    V1 = V(:,1:number_frequecies);  % truncate
    
    Sigma_s = Sigma(1:number_frequecies,1:number_frequecies); % truncate
    H1_hat = U1*Sigma_s*(V1');
    Psi_hat = pinv(H1_hat)*H2;
    eigenvalues = eigs(Psi_hat,number_frequecies); % sorted and truncated to "s" largest values
    frequencies = real(log(eigenvalues)/(-2*pi*1i));

    % transform to positive values and sort due to size 
    frequencies = frequencies + (frequencies<0);
    frequencies= sort(frequencies);
    
end