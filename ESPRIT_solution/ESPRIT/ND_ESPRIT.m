function [frequencies_doppler,frequencies_angle] = ND_ESPRIT(signal,signal_length_doppler,signal_length_angle,number_frequecies)
    
    M_doppler = signal_length_doppler;
    M_angle = signal_length_angle;
    
    
    L_doppler = floor((M_doppler/3+1)); % optimal separation condition 
    L_angle = floor((M_angle/3+1)); % optimal separation condition 
    
    L1 = L_doppler;
    L2= L_angle;
    K_doppler = M_doppler - L_doppler + 1; 
    K_angle = M_angle - L_angle + 1;
    K1=K_doppler;
    K2= K_angle;

    % y is tmp tool to create Hankell tensor
    y_tmp = zeros(K1,K2,L1*L2);
    for i1 = 1:K1
        for i2 = 1:K2         
         vectorized_sub = reshape(signal(i1:i1+L1-1,i2:i2+L2-1),1,[]);
         y_tmp(i1,i2,:) = vectorized_sub';  % (4x1)
        end
    end
    
    % create Hankel tensor
    H = zeros(L1*L2,K1*K2);
    for i1=1:K1
        for i2 =1:K2
            H(:,i2+(i1-1)*K2) = y_tmp(i1,i2,:); 
        end
    end
    
    %SVD 
    [U,Sigma,V] = svd(H);
    U_tilde = U(:,1:number_frequecies);  % truncate
    
    % kronecker product 
    I_L1 = eye(L1);
    I_L2= eye(L2);
 
    I_L1_up = I_L1(2:L1,:);
    I_L1_down = I_L1(1:L1-1,:);
    
    I_L2_up = I_L2(2:L2,:);
    I_L2_down = I_L2(1:L2-1,:);
   
    U_shift_1_up = kron(I_L1_up,I_L2)*U_tilde;
    U_shift_1_down = kron(I_L1_down,I_L2)*U_tilde;
    
    U_shift_2_up = kron(I_L1,I_L2_up)*U_tilde;
    U_shift_2_down = kron(I_L1,I_L2_down)*U_tilde;
    
    % compute F
    F1 = pinv(U_shift_1_down)*U_shift_1_up;
    F2 = pinv(U_shift_2_down)*U_shift_2_up;
    
    % compute K, beta = 1
    K = F1 + F2 ; 
    
    % SVD nr 2
    %[T,Eta,T_inv] = svd(K);
    [T,~] = eig(K);
    
    D1 = T\F1*T;
    D2 = T\F2*T;
    
    modes_doppler = diag(D1);
    modes_angle = diag(D2);
     
    frequencies_doppler = real(log(modes_doppler)/(-2*pi*1i));
    frequencies_angle = real(log(modes_angle)/(-2*pi*1i));
  
    frequencies_doppler = frequencies_doppler + (frequencies_doppler<0);
    frequencies_angle = frequencies_angle + 1*(frequencies_angle<0);  
end