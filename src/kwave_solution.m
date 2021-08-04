%function [p, exec_time] = kwave_solution(L, sos_map, absorption_coeff_map, source_location, omega, min_sos, ~, cfl, roundtrips, dx)
function p = kwave_solution(L, sos_map, absorption_coeff_map, source_location, omega, min_sos, ~, cfl, roundtrips, dx)
    % Solving the 2D Helmholtz equation as a steady state of
    % the wave propagation using kwave
    
    % set the literals
    Nx = double(L);           % number of grid points in the x (row) direction
    Ny = double(L);           % number of grid points in the y (column) direction
    dx = double(dx);
    dy = dx;                  % grid point spacing in the y direction [m]
    w0 = double(omega);
    source_freq = w0/(2*pi);   % [Hz]
    source_mag = 1.; 
    min_sos = double(min_sos);
    cfl = double(cfl);
    roundtrips = double(roundtrips);
    
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
    % define the properties of the propagation medium
    medium.sound_speed = ones(L,L).*sos_map;	% [m/s]
    medium.sound_speed_ref = min(sos_map(:));
    medium.density = 1.;
    medium.alpha_coeff = absorption_coeff_map;
    medium.alpha_power = 0;
    
    % calculate the time step using an integer number of points per period
    ppw = (min_sos/source_freq)/dx;                    % points per wavelength
    ppp = ceil(ppw / cfl);      % points per period
    T   = 1 / source_freq;             % period [s]
    dt  = T / ppp;              % time step [s]
    
    % calculate the number of time steps to reach steady state
    t_end = roundtrips*sqrt( kgrid.x_size.^2 + kgrid.y_size.^2 )/min(medium.sound_speed(:)); % must be 60*
    Nt = round(t_end / dt);
    
    % create the time array
    kgrid.setTime(Nt, dt);
    
    source.p_mask = zeros(L);
    source.p_mask(source_location(1),source_location(2)) = 1;
    
    % define the input signal
    source.p = createCWSignals(kgrid.t_array, source_freq, source_mag, 0);
    
    % set the sensor mask to cover the entire grid
    
    sensor.mask = ones(Nx, Ny);
    % record the last 3 cycles in steady state
    num_periods = 3;
    T_points = round(num_periods * T / kgrid.dt);
    sensor.record_start_index = Nt - T_points + 1;
    
    % input arguments
%    input_args = {'CartInterp', 'nearest','PMLSize',15,'PMLAlpha', 4, ...
%        'DataPath','/tmp/','DataName','sample','DeleteData', false};
    input_args = {'CartInterp', 'nearest','PMLSize',15,'PMLAlpha', 4, 'PlotSim', false};
    

    % run the simulation
    %sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    [amp, phase] = extractAmpPhase(sensor_data, 1/kgrid.dt, source_freq, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    
    p = amp.*exp(1i*phase);
    p = reshape(p, [L,L]);
    
    % Get exec time
%    exec_time = h5readatt(['/tmp/', 'sample_output.h5'], '/', 'simulation_phase_execution_time');
%    exec_time = exec_time(~isspace(exec_time));
%    exec_time(exec_time == 's') = [];
%    exec_time = str2num(exec_time);
    
    end
    