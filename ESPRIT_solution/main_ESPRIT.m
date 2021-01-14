%% INCLUDE PATHS
addpath(genpath('./ESPRIT'))

%% SET CONSTANS
number_targets = 1;
nbr_samples = 40; 
nbr_antenna_elements = 40;

if (nbr_samples ~= nbr_antenna_elements)
    disp('ESPRIT requires the dimension of samples and antenna elements ot be equal!')
end
%% IMPORT DATA
dataStruct = load('filename');
signal = dataStruct.S;

%% Perform ND-ESPRIT
[frequencies_angle,frequencies_doppler] = ND_ESPRIT(signal,nbr_samples,nbr_antenna_elements,number_targets);
 