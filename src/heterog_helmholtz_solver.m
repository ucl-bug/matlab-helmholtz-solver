% Solves the absorbing Helmholtz equation using GMRES

function p = heterog_helmholtz_solver(omega, c0, abs_coeff, dx, N, src_location, pml_size, sigma_star)
    arguments
        omega (1,1) double {mustBeNumeric, mustBeScalar(omega)}
        c0 (:,:) double {mustBeNumeric}
        abs_coeff (:,:) double {mustBeNumeric}
        dx (1,1) double {mustBeNumeric, mustBeScalar(dx)}
        N (1,1) {mustBeNumeric, mustBeScalar(N)}
        src_location (:,1) {mustBeNumeric}
        pml_size (1,1) {mustBeNumeric, mustBeScalar(pml_size)} = 15
        sigma_star (1,1) {mustBeNumeric, mustBeScalar(sigma_star)} = 4*omega
    end 
    
    max_iter = 1000; 

    % set up the source
    source = zeros(N,N);
    source(src_location(1),src_location(2)) = 1/(dx.^2);

    % ensure the sound speed has the right dimensions
    c0 = ones(size(source))*c0;
    
    % Convert absorption coefficient from abs_coeff in dB/cm (as used in k-Wave) to the units required in the Helmholtz solver
    absorption = 100/(20*log10(exp(1)))*abs_coeff.*c0/omega;
    
    % GMRES solution    
    [p, ~] = helmholtz_solver(c0, source, absorption, ...
        omega, pml_size, sigma_star, max_iter, 50, ...
        1e-6, source*0, dx);
end

