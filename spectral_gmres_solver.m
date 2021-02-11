function [p, relres, A, B, M] = spectral_gmres_solver(sos_map, absorption, source, omega, pml_size, ...
        sigma_star, max_iter, checkpoint_frequency, tol, guess, A, B);

    % Spectral Helmholtz solution with attenuation using a PML
    %
    % author: Bradley Treeby, Antonio Stanziola
    % date: 9th April 2020
    % last update: 11th Feb 2021 (Antonio Stanziola)

    % Preparing settings
    assert(...
        size(sos_map, 1) == size(sos_map, 2), ...
        'The speed of sound map must have the same width and height');
    L = size(sos_map, 1);
    Nx = double(L); % number of grid points in the x (row) direction
    Ny = double(L); % number of grid points in the y (column) direction
    dx = double(1); % grid point spacing in the y direction [m]
    dy = dx;
    w0 = double(omega);

    % Builds the Laplacian matrices, if not provided
    if A == 0
        disp("--> Building laplacian operator")
        % Create wavenumbers
        kx = (-pi / dx):2 * pi / (dx * Nx):(pi / dx - 2 * pi / (dx * Nx));
        ky = (-pi / dy):2 * pi / (dy * Ny):(pi / dy - 2 * pi / (dy * Ny));
        kx = ifftshift(reshape(kx, [], 1));
        ky = ifftshift(reshape(ky, 1, []));
        k_sq = kx.^2 + ky.^2;

        % Expand wavenumbers
        kx_mat = repmat(kx, 1, Ny);
        ky_mat = repmat(ky, Nx, 1);

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

        % form 1D DFT matrices
        dftx = fft (eye(Nx));
        dfty = fft (eye(Ny));
        idftx = ifft(eye(Nx));
        idfty = ifft(eye(Ny));

        % Transform everything in tensor ops
        dftx = tensor_ops(dftx, 'matrix_to_LHS_product');
        idftx = tensor_ops(idftx, 'matrix_to_LHS_product');
        dfty = tensor_ops(dfty, 'matrix_to_RHS_product');
        idfty = tensor_ops(idfty, 'matrix_to_RHS_product');
        inv_gamma_x = tensor_ops(1 ./ gamma_x, 'matrix_to_elementwise');
        inv_gamma_y = tensor_ops(1 ./ gamma_y, 'matrix_to_elementwise');
        ikx = tensor_ops(1i * kx_mat, 'matrix_to_elementwise');
        iky = tensor_ops(1i * ky_mat, 'matrix_to_elementwise');

        % Making derivative operators
        Dx = idftx * (ikx * dftx);
        Dy = idfty * (iky * dfty);

        % form A matrix  1/gx d/dx (1/gx d/dx)
        A = inv_gamma_x * (Dx * (inv_gamma_x * Dx));

        % form B matrix  1/gy d/dy (1/gy d/dy)
        B = inv_gamma_y * (Dy * (inv_gamma_y * Dy));
    end

    % form D term
    D = tensor_ops((1 + 1i * absorption) * w0^2 ./ (sos_map.^2), 'matrix_to_elementwise');

    % Full system matrix
    M = A + B + D;

    % Making source
    src_vec = tensor_ops(source, 'matrix_to_vector');

    % Solve problem with GMRES
    disp("--> Solving with GMRES")
    x0 = tensor_ops(guess, 'matrix_to_vector'); %0 .* rand(length(src_vec), 1);
    maxit = checkpoint_frequency;
    restarts = floor(max_iter / checkpoint_frequency);

    [x, flag, relres] = gmres(M, src_vec, restarts, tol, maxit, [], [], x0);
    p = tensor_ops(x, 'vector_to_matrix', [Nx, Ny]);