classdef field2simulation
    properties
       fs = 20e6; %  Sampling frequency [Hz]
       c = 1540;                     %  Speed of sound [m/s]
       width = 1e-3;                %  Width of element [m]
       element_height = 1e-3;        %  Height of element [m]
       kerf = 0.5e-3;                    %  kerf [m]
       focal_point
       n_elements
       n_cycles
       tx
       p_max
       hp_det
       hp_det_norm_dB
       grid
       start_time_h_tx
       h_tx
       start_time_hp
    end
    methods
        function obj = field2simulation(f0,w_array)
            % Initialize Field
            addpath('C:\Users\Admin\Documents\Studies\Masters\Research\Code\Field_II_ver_3_22_windows');
            field_init(0);

            % Create the array transducer
            if nargin < 2
                w_array = 4e-2;             %  Array width [m]
            end
            if nargin < 1
                f0=5e6;                     %  Transducer center frequency [Hz]
            end
            if nargin < 3
                obj.focal_point = [0 0 40]/1000;                                    %  Fixed focal point [m]
            end
            obj.n_elements=floor(w_array/(obj.width+obj.kerf));                 %  Number of physical elements in the transmit aperture
            n_cycles = 2;                                           %  Number of cycles for transmission

            %  Set the relevent simulation parameters
            set_sampling(obj.fs);                   %  Sets sampling frequency
            set_field('use_triangles',0);       %  Tells whether to use triangles (1) or not (0)
            set_field('use_rectangles',1);      %  Tells whether to use rectangles (1) or not (0)
            set_field('use_att',1);             %  Tells whether to use attenuation (1) or not (0)
            set_field('att',80);               % base atteuation [db/m]
            set_field('freq_att', 6e-4);          % frequency dependent atteuation [db/m/Hz]
            set_field('att_f0', 0.5e6);         % minimum attenuation frequency
            set_field('c',obj.c);                   %  Sets the speed of sound

            %  Generate aperture for transmission
            obj.tx = xdc_linear_array(...
                obj.n_elements,...
                obj.width,...
                obj.element_height,...
                obj.kerf,...
                1,...
                1,...
                obj.focal_point);

            %  Set the impulse response of the transmit aperture
            impulse_response=sin(2*pi*f0*(0:1/obj.fs:n_cycles/f0));
            impulse_response=impulse_response.*hanning(length(impulse_response))';

            xdc_impulse(obj.tx,impulse_response);

            %  Set the excitation of the transmit aperture
            excitation=sin(2*pi*f0*(0:1/obj.fs:n_cycles/f0));
            xdc_excitation(obj.tx,excitation);
            
            %% Plot the signal
%             figure;
%             subplot(2,1,1);
%             plot(0:1/obj.fs:2/f0,excitation);
%             title('Excitation Signal');xlabel('t[sec]');
%             subplot(2,1,2);
%             plot(0:1/obj.fs:2/f0,impulse_response);
%             title('Impulse Response');xlabel('t[sec]');

            %%
            % Apodize the array on transmit and/or receive
            xdc_apodization(obj.tx,0,ones(1,obj.n_elements));


            % Define matrix of [x y z] coordinates for where the field calculations occur
            W = 1e-2;                   % [m]
            [x,y,z] = meshgrid(-W:0.0001:W,...
                0,...
                [-2*W:0.0001:2*W]+obj.focal_point(3));
            obj.grid = [x(:), y(:), z(:)];


            % Calculate the emitted field
            [hp, obj.start_time_hp]=calc_hp(obj.tx, obj.grid);
            obj.p_max = reshape(sqrt(sum(hp.^2,1)),...
                [size(x,2),size(x,3)]);             % Max over the y-axis
        end

        function plot_pressure_field(obj)
            %
            x = unique(obj.grid(:,1));
            z = unique(obj.grid(:,3));
            figure;

            % Display the pressure field
            lateral_dim=1000*x; %[mm]
            axial_dim=1000*z;
            P1=rot90(obj.p_max,1);
            Result=flipud(P1);
            imagesc(lateral_dim,axial_dim,dbscale(Result,50));
            colormap(gray);
            xlabel('Lateral distance (mm)');
            ylabel('Depth (mm)');
            title('Transmitted field');
            colorbar;
            %

            % Plot the axial cross-section at the centerline
            figure;
            center_line = floor(size(Result,2)/2)+1;
            subplot(1,2,1);
            plot(axial_dim,Result(:,center_line));
            xlabel('Depth (mm)');
            ylabel('Normalized Pressure');
            title('Axial cross-section');
            [max_axial,I_max_axial] = max(Result(:,center_line));
            fwhm_axial = find_fwhm(Result(:,center_line),...
                find(z == obj.focal_point(3)));
            hold on;
%             xline(axial_dim(...
%                 max([I_max_axial-0.5*fwhm_axial,1])),'r--');
%             xline(axial_dim(...
%                 min([I_max_axial+0.5*fwhm_axial,length(axial_dim)])),'r--');
            hold off;
            
            % Plot the lateral cross-section at the focal distance
            [~, ind] = min(abs(z - obj.focal_point(3)));
            obj.hp_det_norm_dB = convert2db(Result);
            subplot(1,2,2);
            plot(lateral_dim, obj.hp_det_norm_dB(ind,:));
            fwhm_lateral = find_fwhm(Result(ind,:));
            xlabel('Lateral distance (mm)');
            ylabel('Normalized Pressure (dB)');
            title('Lateral cross-section');
            hold on;
            [v_max_lateral,I_max_lateral] = max(Result(ind,:));
%             xline(max([0,I_max_lateral-0.5*fwhm_lateral]),'r--');
%             xline(min([size(Result,2),I_max_lateral+0.5*fwhm_lateral]),'r--');
            hold off;
            
            fprintf('FWHM in axial direction: %f [mm]\n',...
                fwhm_axial/10);

            fprintf('FWHM in lateral direction: %f [mm]\n',...
                fwhm_lateral/10);
        
        end
    end
end