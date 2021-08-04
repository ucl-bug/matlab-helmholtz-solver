function [p, relres] = helmholtz_solver(sos_map, source, absorption, ...
    omega, pml_size, sigma_star, max_iter, checkpoint_frequency, tol, guess, dx)
    % Spectral Helmholtz solver with attenuation using a PML
    %
    % author: Bradley Treeby, Antonio Stanziola, Santeri Kaupinmäk
    % date: 9th April 2020
    % last update: 25th Feb 2021 (Antonio Stanziola)

    % TODO:
    % - Support arbitrary dimensions, now is just 2D
    % - Support anisotropic grids
    
    arguments
        sos_map (:,:) {mustBeNumeric}
        source (:,:) {mustBeNumeric,mustBeEqualSize(source, sos_map)} 
        absorption {mustBeNumeric} = 0.0
        omega (1,1) double {mustBeNumeric, mustBeScalar(omega)} = 1.0
        pml_size (1,1) int64 {mustBeNumeric, mustBeScalar(pml_size)} = 12
        sigma_star (1,1) double {mustBeNumeric, mustBeScalar(sigma_star)} = 2.0
        max_iter (1,1) int64 {mustBeNumeric, mustBeScalar(max_iter)} = 1000
        checkpoint_frequency (1,1) int64 {mustBeNumeric, mustBeScalar(checkpoint_frequency)} = 20
        tol (1,1) double {mustBeNumeric, mustBeScalar(tol)} = 1e-6
        guess (:,:) {mustBeNumeric,mustBeEqualSize(guess, source)} = zeros('like', source)
        dx (1,1) double {mustBeNumeric, mustBeScalar(dx)} = 1.0
    end

    % Temporary validations
    assert(...
        size(sos_map, 1) == size(sos_map, 2), ...
        'The speed of sound map must have the same width and height');

    % Local variables
    L = size(sos_map, 1);
    Nx = double(L); % number of grid points in the x (row) direction
    Ny = double(L); % number of grid points in the y (column) direction
    pml_size = double(pml_size);
    dy = dx;        % grid point spacing in the y direction [m]
    w0 = double(omega);

    %% Derivatives with PML

    % Create wavenumbers
    kx = (-pi / dx):2 * pi / (dx * Nx):(pi / dx - 2 * pi / (dx * Nx));
    ky = (-pi / dy):2 * pi / (dy * Ny):(pi / dy - 2 * pi / (dy * Ny));
    kx = ifftshift(reshape(kx, [], 1));
    ky = ifftshift(reshape(ky, 1, []));

    % Expand wavenumbers
    kx_mat = 1i *repmat(kx, 1, Ny);
    ky_mat = 1i *repmat(ky, Nx, 1);

    % define sigma
    sigma_x = sigma_star .* (((1:pml_size) - pml_size) ./ pml_size).^2;
    sigma_y = sigma_star .* (((1:pml_size) - pml_size) ./ pml_size).^2;
 
    % define PML gamma
    gamma_x = ones(Nx, 1);
    gamma_x(1:pml_size) = 1 + 1i * sigma_x / w0;
    gamma_x(end - pml_size + 1:end) = flip(gamma_x(1:pml_size));
 
    gamma_y = ones(Nx, 1);
    gamma_y(1:pml_size) = 1 + 1i * sigma_y / w0;
    gamma_y(end - pml_size + 1:end) = flip(gamma_y(1:pml_size));
 
    % expand the PML function
    gamma_x = repmat(gamma_x, 1, Ny);
    gamma_y = repmat(transpose(gamma_y), Nx, 1);
    inv_gamma_x = 1./gamma_x;
    inv_gamma_y = 1./gamma_y;

    % Make derivative functions
    Dx = @(field) derivative_with_pml(field, kx_mat, inv_gamma_x, 1);
    Dy = @(field) derivative_with_pml(field, ky_mat, inv_gamma_y, 2);

    %% Vectorize 
    source = source(:);
    guess = guess(:);

    % Helmholtz operator function handler
    helmholtz_op = @(field) helmholtz_pde(field, Dx, Dy, absorption, w0, sos_map, Nx, Ny);

    %% Solve with GMRES
    [p, ~, relres] = gmres(helmholtz_op, source, checkpoint_frequency, tol, max_iter, [], [], guess);
    
    % Make 2D
    p = reshape(p,[Nx Ny]);
end

function y = helmholtz_pde(field, Dx, Dy, absorption, w0, sos_map, Nx, Ny)   
    % Make 2D
    field = reshape(field,[Nx Ny]);
    
    % Make laplacian
    L = Dx(Dx(field)) + Dy(Dy(field));

    % zero-th order therm
    z = ((1 + 1i * absorption) .* w0 ./ (sos_map)).^2 .* field;

    % Return as a vector
    y = L + z;
    y = y(:);
end

function df = derivative_with_pml(field, kvector, inv_gamma, axis)
    arguments
        field (:,:) double 
        kvector (:,:) double {mustBeEqualSize(field, kvector)}
        inv_gamma (:,:) double {mustBeEqualSize(field, inv_gamma)}
        axis (1,1) int64 {mustBeScalar(axis)}
    end

    df = inv_gamma .* spectral_derivative(field, kvector, axis);
end

function df = spectral_derivative(field, kvector, axis)
    arguments
        field (:,:) double 
        kvector (:,:) double {mustBeEqualSize(field, kvector)}
        axis (1,1) int64 {mustBeScalar(axis)}
    end

    if axis == 1
        df = ifft(kvector.*fft(field));
    elseif axis == 2
        % Transpose, derive, transpose again.
        % Use `conj` because matlab doesnt differenciate between
        % transpse and hermitian transpose
        df = conj(spectral_derivative(conj(field)', conj(kvector)', 1))';
    else
        ME = MException('derivative:wrongAxis', ...
        'axis must be 1 or 2, found %s',str(axis));
        throw(ME);
    end
end