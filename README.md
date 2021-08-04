# Matlab Helmholtz Solver
A pedagogical solver for the 2D Helmholtz equation, using a Fourier Pseudospectral method.

## Example
```matlab
%% Setup
N = 128;                                     % number of grid points in each direction
size_x = 10/1000;                            % Domain size [m]
dx = size_x / N;                             % step size [m]

omega = 2*(10^6)*2*pi;                       % angular frequency [rads/s], 2*pi*f_0
c0 = 1300;                                   % Homogeneous speed of sound [m/s]
abs_coeff = zeros(N,N);                      % absorption coefficient [dB/cm]
abs_coeff(N/4:3*N/4,N/4:3*N/4)

src_location = [25,35];                      % Source location in gridpoints
pml_size = 15;                               % PML size in gridpoints
sigma_star = 4*omega;                        % Maximum PML amplitude

% run the GMRES Helmholtz solver
p = heterog_helmholtz_solver(omega, c0, abs_coeff, dx, N, src_location, pml_size, sigma_star);
```

## Documentation
Please visit the [Documentation.md](Documentation.md) file for information on how to use the solver.


## Contributors

[Antonio Stanziola](https://bug.medphys.ucl.ac.uk/antonio-stanziola) - [@astanziola](https://github.com/orgs/ucl-bug/people/astanziola)

[Bradley Treeby](https://bug.medphys.ucl.ac.uk/bradley-treeby) - [@btreeby](https://github.com/orgs/ucl-bug/people/btreeby)

[Santeri Kaupinmäki](https://bug.medphys.ucl.ac.uk/santeri-kaupinmaki)

[Ben Cox](https://bug.medphys.ucl.ac.uk/ben-cox) - [@bencox](https://github.com/orgs/ucl-bug/people/bencox)