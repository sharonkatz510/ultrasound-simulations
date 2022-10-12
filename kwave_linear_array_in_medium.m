classdef kwave_linear_array_in_medium
    properties
        f0
        transducer
        kgrid
        medium
        sensor_data
        p_max
        c = 1540;                           % speed-of-sound [m]
        total_transducer_width              % [m]
        L_element = 15e-3;                  % [m]
        K_element = 0;                      % [m]
        W_element = 0.5e-3;                 % [m]
    end
    methods
        function obj = kwave_linear_array_in_medium(f0, focus, varargin)
            % Constructor for kwave_linear_array_in_medium simulation
            % object
            % Inputs:
            %   - f0 = Transmission frequency [Hz]
            %   - focus = [depth, y] coordinates of the focus [m]
            % Optional inputs:
            %   - n_cycles
            %   - total_transducer_width [m]
            %   - grid = kWave grid object (may be created using the
            %       'kWaveGrid' function
            %   - medium = struct containing the following fields:
            %       sound_speed [m/s]
            %       density [kg/m3
            %       BonA (nonlinearity parameter)
            %       alpha_coeff (power law absorption prefactor)
            %           [dB/(MHz*cm)]
            %       alpha_power (power law absorption exponent)
            addpath(genpath('C:\Program Files\Polyspace\R2021a\toolbox\k-Wave'));
            
            % Check inputs
            P = inputParser();
            P.addParameter('n_cycles',2,...
                @(x) isa(x,'double') && sum(size(x) == [1,1])==2);
            P.addParameter('total_transducer_width',4e-2,@(x) isnumeric(x));
            P.addParameter('grid',[],@(x) isa(x, 'kWaveGrid'));
            P.addParameter('medium',[],@(x) isstruct(x));
            P.parse(varargin{:});
            assert((length(focus)==2) && (size(focus,1)*size(focus,2)==2)...
                && isnumeric(focus),...
                'Expected focus to be a numeric 1x2 vetor');
            n_cycles = P.Results.n_cycles;
            obj.total_transducer_width = P.Results.total_transducer_width;
            obj.f0 = f0;
            
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
            x = 1.25*focus(1);                  % axial axis length [m]
            y = 1.25*(obj.total_transducer_width + focus(2));                   % lateral axis width [m]
            
            % calculate the spacing between the grid points
            dx = min([(obj.c)/(2*f0),... % set to comply with Nyquist theorem for 5 MHz[m]
                obj.W_element]);        % make sure not to exceed the size of trnsducer elemetns         
            dy = dx;                    % [m]
            dz = dx;                    % [m]
            
            % set total number of grid points not including the PML
            Nx = 2^nextpow2(ceil(x/dx));    % [grid points]
            Ny = 2^nextpow2(ceil(y/dy));    % [grid points]
            Nz = 4;     % [grid points]
            
            % create the k-space grid
            if isempty(P.Results.grid)
                obj.kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
            else
                obj.kgrid = P.Results.grid;
            end

            % =========================================================================
            % DEFINE THE MEDIUM PARAMETERS
            % =========================================================================
            if ~isempty(P.Results.medium)
                unique_dims = [];
                required_fields = {'sound_speed','density','BonA',...
                    'alpha_coeff','alpha_power'};
                for field = required_fields
                    assert(sum(contains(fieldnames(P.Results.medium),field))==1,...
                        'Missing medium property');
                    dim = size(getfield(P.Results.medium,field{:}));
                    if dim ~= [1,1]
                        unique_dims = unique([unique_dims, dim]);
                        assert(size(unique_dims,1)<=1,...
                            "medium coefficient dimensions don't match");
                        assert(...
                            sum(dim == [obj.kgrid.Nx, obj.kgrid.Ny, obj.kgrid.Nz])==3,...
                            "medium coefficient dimensions don't match the k-grid");
                    end
                end
                obj.medium = P.Results.medium;
            else                   
                % define the properties of the propagation medium
                obj.medium.sound_speed = 1540;      % [m/s]
                obj.medium.density = 1000;          % [kg/m^3]
                obj.medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
                obj.medium.alpha_power = 1.5;
                obj.medium.BonA = 6;
            end
            
            % =========================================================================
            % DEFINE THE ULTRASOUND TRANSDUCER
            % =========================================================================

            % physical properties of the transducer
            transducer.element_width = max([1,floor(obj.W_element/dy)]);       % width of each element [grid points]
            transducer.element_spacing = floor(obj.K_element/dy);     % spacing (kerf width) between the elements [grid points]
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
            transducer.focus_distance = focus(1);              % focus distance [m]
            transducer.elevation_focus_distance = transducer.position(3)*dz;    % focus distance in the elevation plane [m]
            transducer.steering_angle = rad2deg(atan(focus(2)/focus(1)));                  % steering angle [degrees]

            % apodization
            transducer.transmit_apodization = 'Rectangular';    
            transducer.receive_apodization = 'Rectangular';

            % define the transducer elements that are currently active
            transducer.active_elements = ones(transducer.number_elements, 1);
            
            % =========================================================================
            % DEFINE THE INPUT SIGNAL
            % =========================================================================
            % create the time array
            t_end = (1.5*x/...
                (cos(deg2rad(transducer.steering_angle))*obj.medium.sound_speed))...
                + (n_cycles/f0);                  % [s]
            obj.kgrid.makeTime(obj.medium.sound_speed, [], t_end);
            
            % define properties of the input signal
            source_strength = 10e6;          % [Pa]

            % create the input signal using toneBurst 
            input_signal = toneBurst(1/obj.kgrid.dt, f0, n_cycles);

            % scale the source magnitude by the source_strength divided by the
            % impedance (the source is assigned to the particle velocity)
            input_signal = (source_strength ./...
                (obj.medium.sound_speed * obj.medium.density)) .* input_signal;


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
            [obj.sensor_data] = kspaceFirstOrder3D(obj.kgrid, obj.medium, obj.transducer, sensor, input_args{:});
            obj.p_max = reshape(obj.sensor_data.p_max,Nx,Ny,Nz);
        end
        function plot_pressure_field(obj)
            % plot the final wave-field
            figure;
            imagesc(obj.kgrid.y_vec * 1e3,...
                (obj.kgrid.x_vec- min(obj.kgrid.x_vec)) * 1e3,...
                max(obj.p_max,[],3));
            colormap(getColorMap);
            ylabel('x-position [mm]');
            xlabel('y-position [mm]');
            axis image;
            title('Maximal pressure');


            figure;
            subplot(1,2,1);
            plot(obj.kgrid.y_vec * 1e3,...
                max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],3));
            title('Max pressure at the focal plane (x = 40mm)');
            ylabel('Pressure [Pa]');
            xlabel('y-position [mm]');

            subplot(1,2,2);
            plot(obj.kgrid.x_vec * 1e3,...
                max(obj.p_max(:,round(obj.kgrid.Ny/2),:),[],3));
            title('Max pressure at the center-line');
            ylabel('Pressure [Pa]');
            xlabel('x-position [mm]');
            
            fprintf('FWHM in lateral direction: %f [mm]\n',...
                find_fwhm(...
                max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],3))*...
                obj.kgrid.dy*1e3);

            fprintf('FWHM in axial direction: %f [mm]\n',...
                find_fwhm(...
                max(obj.p_max(:,round(obj.kgrid.Ny/2),:),[],3))...
                *obj.kgrid.dx*1e3);
            
            fprintf('Maximal pressure: %f [precent of transmission]\n',...
                100*max(obj.p_max(round(40e-3/obj.kgrid.dx),:,:),[],'all')./...
                ((10e6)*obj.transducer.number_elements));
        end
    end
end