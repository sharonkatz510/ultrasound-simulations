classdef kwave_linear_array_in_medium
    properties
        transducer
        kgrid
        sensor_data
        p_max
        c = 1540;
        total_transducer_width = 4e-2;      % [m]
        L_element = 1e-3;                   % [m]
        K_element = 0.5e-3;                 % [m]
        W_element = 1e-3;                   % [m]
    end
    methods
        function obj = kwave_linear_array_in_medium(f0, n_cycles)
            addpath(genpath('C:\Program Files\Polyspace\R2021a\toolbox\k-Wave'));
            % simulation settings
            DATA_CAST = 'single';

            % =========================================================================
            % DEFINE THE K-WAVE GRID
            % =========================================================================

            % set the size of the perfectly matched layer (PML)
            PML_X_SIZE = 20;            % [grid points]
            PML_Y_SIZE = 10;            % [grid points]
            PML_Z_SIZE = 10;            % [grid points]

            % set desired grid size in the x-direction not including the PML
            x = 50e-3;                  % axial axis length [m]
            y = 4e-2;                   % lateral axis width [m]
            
            % calculate the spacing between the grid points
            dx = min([(obj.c)/(2*f0),... % set to comply with Nyquist theorem [m]
                obj.W_element]);        % make sure not to exceed the size of trnsducer elemetns         
            dy = dx;                    % [m]
            dz = dx;                    % [m]
            
            % set total number of grid points not including the PML
            Nx = 2^nextpow2(ceil(x/dx));    % [grid points]
            Ny = 2^nextpow2(ceil(y/dy));    % [grid points]
            Nz = 4;     % [grid points]
            % create the k-space grid
            obj.kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

            % =========================================================================
            % DEFINE THE MEDIUM PARAMETERS
            % =========================================================================

            % define the properties of the propagation medium
            medium.sound_speed = 1540;      % [m/s]
            medium.density = 1000;          % [kg/m^3]
            medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
            medium.alpha_power = 1.5;
            medium.BonA = 6;

            % create the time array
            t_end = 40e-6;                  % [s]
            obj.kgrid.makeTime(medium.sound_speed, [], t_end);

            % =========================================================================
            % DEFINE THE INPUT SIGNAL
            % =========================================================================

            % define properties of the input signal
            source_strength = 1e6;          % [Pa]

            % create the input signal using toneBurst 
            input_signal = toneBurst(1/obj.kgrid.dt, f0, n_cycles);

            % scale the source magnitude by the source_strength divided by the
            % impedance (the source is assigned to the particle velocity)
            input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;

            % =========================================================================
            % DEFINE THE ULTRASOUND TRANSDUCER
            % =========================================================================

            % physical properties of the transducer
            transducer.element_width = max([1,floor(obj.W_element/dz)]);       % width of each element [grid points]
            transducer.element_spacing = max([1,floor(obj.K_element/dz)]);     % spacing (kerf width) between the elements [grid points]
            transducer.number_elements = ...
                floor(obj.total_transducer_width/...
                (obj.W_element+obj.K_element)) ;    % total number of transducer elements
            
            transducer.element_length = 1;     % length of each element [grid points]
            transducer.radius = inf;            % radius of curvature of the transducer [m]

            % calculate the width of the transducer in grid points
            transducer_width = transducer.number_elements * transducer.element_width ...
                + (transducer.number_elements - 1) * transducer.element_spacing;    % total width in grid points

            % use this to position the transducer in the middle of the computational grid
            transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

            % properties used to derive the beamforming delays
            transducer.sound_speed = 1540;                  % sound speed [m/s]
            transducer.focus_distance = 40e-3;              % focus distance [m]
            transducer.elevation_focus_distance = transducer.position(2)*dx;    % focus distance in the elevation plane [m]
            transducer.steering_angle = 0;                  % steering angle [degrees]

            % apodization
            transducer.transmit_apodization = 'Rectangular';    
            transducer.receive_apodization = 'Rectangular';

            % define the transducer elements that are currently active
            transducer.active_elements = ones(transducer.number_elements, 1);

            % append input signal used to drive the transducer
            transducer.input_signal = input_signal;

            % create the transducer using the defined settings
            obj.transducer = kWaveTransducer(obj.kgrid, transducer);

            % print out transducer properties
            obj.transducer.properties;

            % =========================================================================
            % DEFINE SENSOR MASK
            % =========================================================================

            % create a binary sensor mask with four detection positions
            sensor.mask = ones(Nx,Ny,Nz);
            sensor.record = {'p_final', 'p_max'};

            % =========================================================================
            % RUN THE SIMULATION
            % =========================================================================

            % set the input settings
            input_args = {'DisplayMask', 'off', ...
                'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
                'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};

            % run the simulation
            [obj.sensor_data] = kspaceFirstOrder3D(obj.kgrid, medium, obj.transducer, sensor, input_args{:});
            obj.p_max = reshape(obj.sensor_data.p_max,Nx,Ny,Nz);
        end
        function plot_pressure_field(obj)
            % plot the final wave-field
            figure;
            imagesc(obj.kgrid.y_vec * 1e3,...
                [0:(length(obj.kgrid.x_vec)-1)]* obj.kgrid.dx * 1e3,...
                max(obj.p_max,[],3));
            colormap(getColorMap);
            ylabel('x-position [mm]');
            xlabel('y-position [mm]');
            axis image;
            title('Maximal pressure');


            figure;
            plot(obj.kgrid.y_vec * ...
                1e3,max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],3));
            max_amp = max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],'all');
            hold on;
            yline(max_amp/2);
            legend('Max pressure at the focal plane (x = 40mm)','1/2 P-max');
            ylabel('Pressure [Pa]');
            xlabel('y-position [mm]');

            fprintf('FWHM in lateral direction: %f [mm]\n',...
                find_fwhm(max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],3))...
                *obj.kgrid.dx*1e3);

            fprintf('FWHM in axial direction: %f [mm]\n',...
                find_fwhm(...
                max(...
                squeeze(obj.p_max(round(20e-3/obj.kgrid.dx):end,...
                round(obj.kgrid.Ny/2),:)),[],2)...
                )...
                *obj.kgrid.dy*1e3);
            fprintf('Maximal pressure: %f [Pa]\n', max(obj.p_max,[],'all'));
        end
    end
end