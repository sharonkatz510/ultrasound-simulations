classdef Marmottant_simulation
    properties
        c = 1540;                                 % speed of sound in m/s
        w
        P0
        R0
        rhoL = 1000;                              % liquid density  kg/m^3
        xi = 0.38;                     % elastic compression modulus in N/m
        muL = 0.001;                           % liquid viscosity
        sigma0 = 0.0;                   % initial surface tension in N/m
        sigma_w = 0.074;                 % water surface tension in N/m
        gamma = 1.07;                               % polytropic gas exponent
        kappaS = 2.4e-9 ;                 % shell surface dilatational viscosity in kg/s
        Rbuckling
        Rruptured
        broken
    end
    methods
        function obj = Marmottant_simulation(f0, varargin)
            % Model constructor
            %%%%%%%%%%%%%%%%  constants and parameters     %%%%%%%%%%%%%%%%%%
            P = inputParser();
            P.addOptional('n_cycles',2,...
                @(x) isdouble(x) && sum(size(x)==[1,1])==2);
            P.addOptional('R0',0.75e-6,...
                @(x) isdouble(x) && sum(size(x)==[1,1])==2);
            P.addOptional('P0',1.01e5,...
                @(x) isdouble(x) && sum(size(x)==[1,1])==2);
            P.addOptional('sampling_pressures',[],@(x) isnumeric(x));
            P.parse(varargin{:})
            
            R0 = P.Results.R0;                   % bubble initial radius in meter
            P0 = P.Results.P0;                             % ambient pressure in Pa
            n_cycles = P.Results.n_cycles;                  % number of cycles in waveform

            obj.P0 = P0;
            obj.R0 = R0;
            sampling_pressures = 10.^[1:6];

            obj.Rbuckling = R0;
            obj.Rruptured = 1.1*R0;

            phi = pi;                                 % initial phase of the driving pulse 
            fs = 1e9;                               % sampling frequency in Hz
            %**************   end of parameters     ***********************


            %% ************* Excitation pulse ***************
            obj.w = 2*pi*f0;	                          % angular frequency

            t = 0:1/fs:1.5*n_cycles/f0;
            Ns = length(t);
            sig = sin(2*pi*f0*t+phi).*[zeros(1,100),tukeywin(Ns-200,0.85)',zeros(1,100)];

            expansion = zeros(length(sampling_pressures),length(t));
          
            for n = 1:length(sampling_pressures)

                Pa = sampling_pressures(n);
                filt_pr = Pa*sig; 

                Plist_n = -min(filt_pr);              	% Max negative pressure (as a positive value)
                pdriv = struct('hyd',filt_pr,'t',t,'Pmin',Plist_n);

                obj.broken = 0;
                assignin('base','calling_obj',obj);
                [t_rp,y_rp] = ode45('Marmottant_model',[0 max(pdriv.t)],[R0; 0;],[],pdriv);
                obj = evalin('base','calling_obj');
                r = y_rp(:,1);                   	% wall radius
                rdot = y_rp(:,2);              		% wall velocity

                r_interp = interp1(t_rp,r,t,'linear','extrap');              % interpolate to the same time scale
                rdot_interp = interp1(t_rp,rdot,t,'linear','extrap');        % interpolate to the same time scale
                figure
                plot(t*1e6,r_interp/R0)
                hold on
                plot([t(1) t(end)]*1e6,[obj.Rruptured obj.Rruptured]/R0,'--k')
                xlabel('time (\mus)')
                ylabel('expansion')
                set(gca,'fontsize',14)
                expansion(n,:) = r_interp;

            end

            max_exp = ones(size(sampling_pressures));
            max_exp2 = ones(size(sampling_pressures));
            for n = 1:length(sampling_pressures)
                FF = expansion(n,:);
                max_exp(n) = max(FF/R0);
            end

            figure
            semilogx(sampling_pressures/1e3,max_exp,'r')
            xlabel('Peak negative pressure [kPa]')
            ylabel('Expansion ratio')
            set(gca,'fontsize',14)
        end
    end
end

