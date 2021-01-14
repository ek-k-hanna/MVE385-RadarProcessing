function pointMax = lagrange_1d(S, ri, ci)
% Lagrange over fixed column ---> refine doppler
p_left_doppler = [ri-1, S(ri-1, ci)];   % [ri, ci-1]
p_mid_doppler = [ri, S(ri, ci)];       % [ri,ci]
p_right_doppler = [ri+1, S(ri+1, ci) ] ;  % [ri, ci+1]

% create Lagrange polynomial coefficients
l_left_doppler = ( 1/( (p_left_doppler(1)-p_mid_doppler(1))*(p_left_doppler(1)-p_right_doppler(1)) )) * ...
            [1 -p_mid_doppler(1)-p_right_doppler(1) p_mid_doppler(1)*p_right_doppler(1)];    
l_mid_doppler = ( 1/( (p_mid_doppler(1)-p_left_doppler(1))*(p_mid_doppler(1)-p_right_doppler(1)) )) * ...
            [1 -p_left_doppler(1)-p_right_doppler(1) p_left_doppler(1)*p_right_doppler(1)];   
l_right_doppler = ( 1/( (p_right_doppler(1)-p_mid_doppler(1))*(p_right_doppler(1)-p_left_doppler(1)) )) * ...
            [1 -p_mid_doppler(1)-p_left_doppler(1) p_mid_doppler(1)*p_left_doppler(1)];
lagrange_doppler = p_left_doppler(2)*l_left_doppler + p_mid_doppler(2)*l_mid_doppler + p_right_doppler(2)* l_right_doppler ;
doppler_refined_index = -lagrange_doppler(2)/(2*lagrange_doppler(1));

% Lagrange over fixed row  ---> refine angle
p_left_angle = [ci-1, S(ri, ci-1)];   % [ri, ci-1]
p_mid_angle = [ci, S(ri, ci)];       % [ri,ci]
p_right_angle = [ci+1, S(ri, ci+1) ] ;  % [ri, ci+1]

% create Lagrange polynomial coefficients
l_left_angle = ( 1/( (p_left_angle(1)-p_mid_angle(1))*(p_left_angle(1)-p_right_angle(1)) )) * ...
            [1 -p_mid_angle(1)-p_right_angle(1) p_mid_angle(1)*p_right_angle(1)];     
l_mid_angle = ( 1/( (p_mid_angle(1)-p_left_angle(1))*(p_mid_angle(1)-p_right_angle(1)) )) * ...
            [1 -p_left_angle(1)-p_right_angle(1) p_left_angle(1)*p_right_angle(1)];      
l_right_angle = ( 1/( (p_right_angle(1)-p_mid_angle(1))*(p_right_angle(1)-p_left_angle(1)) )) * ...
            [1 -p_mid_angle(1)-p_left_angle(1) p_mid_angle(1)*p_left_angle(1)];
lagrange_angle = p_left_angle(2)*l_left_angle + p_mid_angle(2)*l_mid_angle + p_right_angle(2)* l_right_angle ;
angle_refined_index = -lagrange_angle(2)/(2*lagrange_angle(1));

pointMax = [doppler_refined_index; angle_refined_index];
end