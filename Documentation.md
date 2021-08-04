# Documentation

## `helmholtz_solver`

Solves the 2D absorbing Helmholtz equation using GMRES. The form of the Helmhotlz equation is

$$
\left[ \nabla^2 + \left(\frac{(1 + i\alpha)\omega}{c_0}\right)^2\right]p(x) = s(x)
$$

No preconditioning is used.

#### Syntax
```matlab
[p, relres] = helmholtz_solver(sos_map, source, absorption, ...
    omega, pml_size, sigma_star, max_iter, checkpoint_frequency, tol, guess, dx)
```

*Parameters*

|  Parameter | Kind | Description |
|---|---| --- |
| `sos_map` | 2D array, scalar | The speed of sound map |
| `source` | 2D array | The source term |
| `absorption` | 2D array, scalar (default=0) | The absorption term |
| `omega` | scalar (default=1) | The angular frequency |
| `pml_size`| scalar (default=12) | The size of the PML in pixels |
| `sigma_star` | scalar (default=2) | The maximum amplitude of the PML $\sigma$ function |
| `max_iter` | scalar (default=1000) | The maximum number of iterations for the GMRES solver |
| `checkpoint_frequency` | scalar (default=1000) | The size of the Krylov subspace between restarts |
| `tol` | scalar (default=1e-6) | The tolerance for the GMRES solver |
| `guess` | 2D array (default=`zeros(N,N)`) | The initial guess for the solution |
| `dx` | scalar (default=1) | The spacing between pixels |

*Outputs*
`p`: The solution to the Helmholtz equation
`relres`: The relative residual of the GMRES solver

## `kwave_solution`

Solves the 2D absorbing Helmholtz equation as the steady state of the wave equation, using [`k-Wave`](http://www.k-wave.org/). Note that [`k-Wave`](http://www.k-wave.org/) must be installed and added to the MATLAB path.

#### Syntax
```matlab
 p = kwave_solution(L, sos_map, absorption_coeff_map, ...
    source_location, omega, min_sos, ~, cfl, roundtrips, dx)
```

*Parameters*

|  Parameter | Kind | Description |
|---|---| --- |
| `L` | 2D array, scalar | The size of the domain |
| `sos_map` | 2D array, scalar | The speed of sound map |
| `absorption_coeff_map` | 2D array, scalar | The absorption coefficient map |
| `source_location` | 1d array of len=2 | The location of the source term |
| `omega` | scalar | The angular frequency |
| `min_sos` | scalar | The minimum speed of sound, for stability calculations |
| `cfl` | scalar  | The CFL number |
| `roundtrips` | scalar | The number of roundtrips to run to consider the wavefield in steady state |
| `dx` | scalar | The spacing between pixels |

*Output*
`p`: The solution to the Helmholtz equation