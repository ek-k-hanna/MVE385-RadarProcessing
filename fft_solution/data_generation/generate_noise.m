function data = generate_noise(amp, samples, antennas)
    % GAUSSIAN NOISE
    gaussian_noise = randn(samples, antennas);
    gaussian_noise = gaussian_noise - mean(gaussian_noise, 'all');
    biased = 1;
    gaussian_noise = gaussian_noise/std(gaussian_noise, biased, 'all');
    
    % UNIFORMAL NOISE WITH REJECTION SAMPLING
    uniformal_noise = rand(samples, antennas);
    
    amp = sqrt(10^(amp/20));
    
    % EULER'S FORMULA
    data = amp*gaussian_noise.*exp(1i*uniformal_noise*pi);
end