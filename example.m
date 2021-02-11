% Setup problem
N = 128;                       % Assumes dx=dy=1
source = zeros(N);
source(32,32) = 1.;
omega = 1;
pml_size = 15;
sigma_star = 4;                % Amplitude of the PML layer at the edges of the domain
max_iter = 1000;               % Max iterations for GMRES
checkpoint_frequency = 20;     % Size of the Krylov subspace for each restart
guess = -source;               % Solution guess
tol = 1e-6;                    % Convergence magnitude for residual in GMRES
A = 0;                         % Set A and B to 0 to automatically calculate the matrix
B = 0;                         %    associated with the Laplacian operator

%{ 
Free space example
%}
disp("Solving homogeneous problem...")
sos_map = ones(N);
attenuation_map = zeros(N);
[p, relres, A, B, M] = spectral_gmres_solver(...
    sos_map, attenuation_map, source, omega, pml_size, sigma_star, ...
    max_iter, checkpoint_frequency, tol, guess, A, B);

% Save image
h = figure('Visible', 'off');
imagesc(real(p));
caxis([-0.05, 0.05]);
colorbar();
print(h,'-dpng','free_space');
close(h);

%{
Heterogeneous example.

Note that the previous example returns the sparse matrices 
A and B that form a decomposition of the laplacian L, such that L = A+B.
They can be passed to the solver to be reused (saves computation)
%}
disp("Solving heterogeneous problem...")
sos_map = ones(N);
sos_map(64:96,64:96) = 1.3;
attenuation_map = zeros(N);
attenuation_map(20:64,64:96) = 0.05;
guess = p;                             % Previous solution should be a good starting point
[p, relres, A, B, M] = spectral_gmres_solver(...
    sos_map, attenuation_map, source, omega, pml_size, sigma_star, ...
    max_iter, checkpoint_frequency, tol, guess, A, B);

% Save image
h = figure('Visible', 'off');
imagesc(real(p));
caxis([-0.05, 0.05]);
colorbar();
print(h,'-dpng','heterogeneous');
close(h);