% Comparison between k-Wave and the GMRES Helmholtz solver in the case of
% heterogeneous absorption


%% Setup
N = 128;                                     % number of grid points in each direction
size_x = 10/1000;                            % Domain size [m]
dx = size_x / N;                             % step size [m]

omega = 2*(10^6)*2*pi;                       % angular frequency [rads/s], 2*pi*f_0
c0 = 1300;                                   % Homogeneous speed of sound [m/s]
abs_coeff = 40*makeDisc(N, N, N/2, N/2, 20); % absorption coefficient [dB/cm]

src_location = [32,32];                      % Source location in gridpoints
pml_size = 15;                               % PML size in gridpoints
sigma_star = 4*omega;                        % chosen for reasons of stability?

% This is run by example_accuracy_test_hetero_abs.m

% run the GMRES Helmholtz solver
p = heterog_helmholtz_solver(omega, c0, abs_coeff, dx, N, src_location, pml_size, sigma_star);

% run the k-Wave solver for a single frequency 
p_kwave = kwave_solution(N, c0, abs_coeff, src_location, omega, c0, 0.0, 0.1, 2, dx);

% Fix amplitude in kwave solution
p_kwave = conj(p_kwave);
p_kwave = p_kwave*p(src_location(1),src_location(2))/ ...
    p_kwave(src_location(1),src_location(2)); 


% PLOTS
figure;
subplot(2,2,1);
imagesc(real(p))
axis off
axis equal
caxis([-0.2, 0.2])
title("GMRES result")
colorbar
axis image

subplot(2,2,2);
imagesc(real(p_kwave))
axis off
axis equal
title("kWave solution")
caxis([-0.2, 0.2])
colorbar
axis image

subplot(2,2,3);
imagesc(100*abs(p - p_kwave)/max(abs(p_kwave(:))))
axis off
axis equal
title("Percent Error")
colorbar
axis image
caxis([0, 2])
drawnow

subplot(2,2,4);
imagesc(abs_coeff)
axis off
axis equal
title("absorption coefficient [dB/cm]")
colorbar
axis image
drawnow