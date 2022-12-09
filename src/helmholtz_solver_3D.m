function [p, relres] = helmholtz_solver_3D(sos_map, source, absorption, ...
    omega, pml_size, sigma_star, max_iter, checkpoint_frequency, tol, guess, dx)
    % Spectral Helmholtz solver with attenuation using a PML
    %
    % author: Bradley Treeby, Antonio Stanziola, Santeri Kaupinmäki
    % date: 9th April 2020
    % last update: 7th Sept 2022 (Santeri Kaupinmaki)

    % TODO:
    % - Support arbitrary dimensions, now is just 3D
    % - Support anisotropic grids
    % - Fix arguments checking
    
%     arguments
%         sos_map (:,:,:) {mustBeNumeric}
%         source (:,:,:) {mustBeNumeric,mustBeEqualSize(source, sos_map)} 
%         absorption {mustBeNumeric} = 0.0
%         omega (1,1) double {mustBeNumeric, mustBeScalar(omega)} = 1.0
%         pml_size (1,1) int64 {mustBeNumeric, mustBeScalar(pml_size)} = 12
%         sigma_star (1,1) double {mustBeNumeric, mustBeScalar(sigma_star)} = 2.0
%         max_iter (1,1) int64 {mustBeNumeric, mustBeScalar(max_iter)} = 1000
%         checkpoint_frequency (1,1) int64 {mustBeNumeric, mustBeScalar(checkpoint_frequency)} = 20
%         tol (1,1) double {mustBeNumeric, mustBeScalar(tol)} = 1e-6
%         guess (:,:,:) {mustBeNumeric,mustBeEqualSize(guess, source)} = zeros('like', source)
%         dx (1,1) double {mustBeNumeric, mustBeScalar(dx)} = 1.0
%     end

    % Temporary validations
    assert(...
        size(sos_map, 1) == size(sos_map, 2), ...
        'The speed of sound map must have the same width and height'); % Ignoring 3rd dimension

    % Local variables
    L = size(sos_map, 1);
    Nx = double(L); % number of grid points in the x (row) direction
    Ny = double(L); % number of grid points in the y (column) direction
    Nz = double(L); % number of grid points in the z (depth) direction
    pml_size = double(pml_size);
    dy = dx;        % grid point spacing in the y direction [m]
    dz = dx;        % grid point spacing in the z direction [m]
    w0 = double(omega);

    %% Derivatives with PML

    % Do in cases where N is even or odd
    if rem(Nx,2)==0
        nx = ((-Nx/2:Nx/2-1)/Nx);
        ny = ((-Ny/2:Ny/2-1)/Ny);
        nz = ((-Nz/2:Nz/2-1)/Nz);
    else
        nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx);
        ny = ((-(Ny-1)/2:(Ny-1)/2)/Ny);
        nz = ((-(Nz-1)/2:(Nz-1)/2)/Nz);
    end
    
    % Create wavenumbers
    kx = ifftshift((2*pi/dx).*nx);
    ky = ifftshift((2*pi/dy).*ny);
    kz = ifftshift((2*pi/dz).*nz);

    [ky_mat, kx_mat, kz_mat] = meshgrid(1i*kx,1i*ky,1i*kz);
    
    % define sigma
    sigma_x = sigma_star .* (((1:pml_size) - pml_size) ./ pml_size).^2;
    sigma_y = sigma_star .* (((1:pml_size) - pml_size) ./ pml_size).^2;
    sigma_z = sigma_y;
 
    % define PML gamma
    gamma_x = ones(Nx, 1);
    gamma_x(1:pml_size) = 1 + 1i * sigma_x / w0;
    gamma_x(end - pml_size + 1:end) = flip(gamma_x(1:pml_size));
 
    gamma_y = ones(Nx, 1);
    gamma_y(1:pml_size) = 1 + 1i * sigma_y / w0;
    gamma_y(end - pml_size + 1:end) = flip(gamma_y(1:pml_size));
    
    gamma_z = ones(Nx, 1);
    gamma_z(1:pml_size) = 1 + 1i * sigma_z / w0;
    gamma_z(end - pml_size + 1:end) = flip(gamma_z(1:pml_size));
 
    % expand the PML function
    [gamma_y, gamma_x, gamma_z] = meshgrid(gamma_x,gamma_y,gamma_z);
    
    inv_gamma_x = 1./gamma_x;
    inv_gamma_y = 1./gamma_y;
    inv_gamma_z = 1./gamma_z;

    % Make derivative functions
    Dx = @(field) derivative_with_pml(field, kx_mat, inv_gamma_x, 1);
    Dy = @(field) derivative_with_pml(field, ky_mat, inv_gamma_y, 2);
    Dz = @(field) derivative_with_pml(field, kz_mat, inv_gamma_z, 3);

    %% Vectorize 
    source = source(:);
    
    guess = guess(:);

    % Helmholtz operator function handler
    helmholtz_op = @(field) helmholtz_pde_3D(field, Dx, Dy, Dz, absorption, w0, sos_map, Nx, Ny, Nz);

    %% Solve with GMRES
    [p, ~, relres] = gmres(helmholtz_op, source, checkpoint_frequency, tol, max_iter, [], [], guess);
    
    % Make 3D
    p = reshape(p,[Nx Ny Nz]);
end

function y = helmholtz_pde_3D(field, Dx, Dy, Dz, absorption, w0, sos_map, Nx, Ny, Nz)   
    % Make 3D
    field = reshape(field,[Nx Ny Nz]);
    
    % Make laplacian
    L = Dx(Dx(field)) + Dy(Dy(field)) + Dz(Dz(field));

    % zero-th order therm
    z = ((1 + 1i * absorption) .* w0 ./ (sos_map)).^2 .* field;

    % Return as a vector
    y = L + z;
    y = y(:);
end

% Derivative with PML and homogeneous density
function df = derivative_with_pml(field, kvector, inv_gamma, axis)
    arguments
        field (:,:,:) double 
        kvector (:,:,:) double {mustBeEqualSize(field, kvector)}
        inv_gamma (:,:,:) double {mustBeEqualSize(field, inv_gamma)}
        axis (1,1) int64 {mustBeScalar(axis)}
    end

    df = inv_gamma .* spectral_derivative(field, kvector, axis);
end

function df = spectral_derivative(field, kvector, axis)
    arguments
        field (:,:,:) double 
        kvector (:,:,:) double {mustBeEqualSize(field, kvector)}
        axis (1,1) int64 {mustBeScalar(axis)}
    end

    if axis == 1
        df = ifft(kvector.*fft(field));
    elseif axis == 2
        df = ifft(kvector.*fft(field, [], 2), [], 2);
    elseif axis == 3
        df = ifft(kvector.*fft(field, [], 3), [], 3);
    else
        ME = MException('derivative:wrongAxis', ...
        'axis must be 1, 2, or 3, found %s',str(axis));
        throw(ME);
    end
end
