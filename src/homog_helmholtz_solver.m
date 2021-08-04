function p = homog_helmholtz_solver(omega, c0, dx, N, src_location, pml_size, sigma_star)
    arguments
        omega (1,1) double {mustBeNumeric, mustBeScalar(omega)}
        c0 (1,1) double {mustBeNumeric, mustBeScalar(c0)}
        dx (1,1) double {mustBeNumeric, mustBeScalar(dx)}
        N (1,1) {mustBeNumeric, mustBeScalar(N)}
        src_location (:,1) {mustBeNumeric}
        pml_size (1,1) {mustBeNumeric, mustBeScalar(pml_size)} = 15
        sigma_star (1,1) {mustBeNumeric, mustBeScalar(sigma_star)} = 4*omega
    end 
    
    max_iter = 1000; 
    
    % GMRES solution
    source = zeros(N,N);
    source(src_location(1),src_location(2)) = 1/(dx.^2);
    c0 = ones(size(source))*c0;
    [p, ~] = helmholtz_solver(c0, source, c0*0.0, ...
        omega, pml_size, sigma_star, max_iter, 50, ...
        1e-6, source*0, dx);
end