function lens_params(sourceloc,dist_mode,num_model,varargin)

global c G pc t_H eps_ % mks universal constants; single precision machine eps
% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
t_H = 3.09E17/0.8; % ! h=0.6 in Schneider
% machine epsilon
eps_ = double(eps('single'));

% unpack varargin
if(nargin == 3)
    flag_plot    = 0;  % plot time delay surf, cross-sections
    flag_summary = 0;  % plot scaled tau(x), mu(x), rho(x) as source moves
    scale_mode   = 1;  % if 1, use scaled lens params; if 0 use mks units
elseif(nargin >= 4)
    flag_plot    = varargin{1};
    flag_summary = varargin{2};
    scale_mode   = varargin{3};
end % if

%% Distances in m between source, lens, observer
switch(dist_mode{1})
    case('z')           % redshifts to source and lens
        [Dd,Ds,Dds] = tdelay(dist_mode);
    case('dist');       % distances in pc
        Dd = dist_mode{2}*pc; Ds = dist_mode{3}*pc;
        Dds = Ds - Dd;  % assume d << D_H
end % switch
sigma_cr = c^2/(4*pi*G) * Ds/(Dd*Dds); % kg/m^2

%% Set up lens plane and convolution
num = pow2(8); % num pixels in [0,max_radius]
max_radius = 25; % box edge
pix = max_radius / num  % pixel size in scaled units
% Need to set R,N so that fns converge to 0 for each SMD:
% kappa(R) = 0 and HT[kappa(R)] = sigma_fft(V) = 0
H = hankel_matrix(0,max_radius,num);
rho = H.r(:); x = [-flipud(H.r(:)); H.r(:)];
if(flag_plot); [x1,x2] = ndgrid(x,x); end; % if plotting plane, make grid

%% Lens parameters
%theta_0 = .009; % ang. radius of Einstein ring in arcsec.
%M = theta_0^2 * (Dd*Ds)/Dds * c^2/(4*G);
M = 1e4*M_earth;                    % total lens mass
mstr = num2str(M/M_sol,'%6.4e');    % string for plotting
r_E = sqrt(4*G*M/c^2 * Dd*Dds/Ds);  % Einstein radius
zeta_0 = r_E;                       % lens plane: x = zeta/zeta_0

% assign parameters for convergence / analytical solns.
switch(num_model)
    case{1}
        model = struct('scale','Einstein radii','title','Schwarzschild lens', ...
            'filestr','swzschild','paramstr',['Mass ',mstr,' M_{\odot}']);
        params = struct('M',M,'r_E',r_E);
        
    case{2}
        if(~scale_mode)
            r_max = r_E;
            sigma_0 = M/(pi*(r_max)^2); % const. density
            x_max = r_max/zeta_0;
            kappa_0 = sigma_0 / sigma_cr;
        else
            x_max = .1             % r_max/zeta_0
            kappa_0 = x_max^(-2)    % scaled density
            % kappa_0 and x_max are NOT independent!
        end % if
        model = struct('scale','Einstein radii','title','Homogeneous disc', ...
            'filestr',strcat('hdisc_x_',stripdec(x_max,'.')), ...
            'paramstr',['Mass ',mstr,' M_{\odot}',' x_0: ',num2str(x_max,'%6.4e')]);
        params = struct('M',M,'r_E',r_E,'x_max',x_max);
        
    case{3} % NFW
        % assume z << 1 so rho_crit equal to present-day value
        rho_crit = 3*t_H^(-2)/(8*pi*G); % ~ 9e-27 kg/m^3 (Hobson)
        % "virial radius" r_200 (nfw-cdm-halos p5)
        r_200 = (M/(200*rho_crit*4*pi/3))^(1/3);
        if(~scale_mode)
            % scale radius r_s, concentration c
            r_s = 10^(-4.5)*pc;     % cutoff radiusof lens
            x_max = r_s/zeta_0      % scaled radius of lens
            c_nfw = r_200/r_s
        else % param. is c
            c_nfw = 2e6             % c > 1 (wright-nfw-halos p3)
            x_max = .5              % turnover radius r_max/r_s
        end
        % characteristic overdensity delta_c (nfw-cdm-halos p5)
        delta_c = (200/3) * c_nfw^3 / (log(1+c_nfw) - c_nfw/(1+c_nfw));
        kappa_s = delta_c*rho_crit*r_200/c/sigma_cr;
        model = struct('scale','Scale radii r_{200}/c', ...
            'title','Navarro-Frenck-White profile', ...
            'filestr',strcat('nfw_x_',stripdec(x_max,'.'),'_c_',stripdec(c_nfw,',')), ...
            'paramstr',['Mass ',mstr,' M_{\odot}',' turnover radius: ',num2str(x_max,'%6.4g'),' c: ',num2str(c_nfw,'%6.4g')]);
        params = struct('M',M,'r_E',r_E,'x_max',x_max,'c',c);
        clear delta_c rho_crit r_s
end % switch

lens.params = params;               % save model parameters
lens.names = model;                 % save model details
clear model params


%% Calculate image locations, magnfication, time delays using HT code
for i = 1:params_len(1)
    params{1} = lens.params.M(i);
    for j = 1:params_len(2)
        params{2} = 
        for k = 1:params_len(3)
            for l = 1:params_len(4)
                params = lens.params{};
                lens = lens_images(sourceloc,num_model,params,varargin);
                % is it necessary to store all data?
                tau(i,j,k,l,:,:) = cellextract(lens.tau_im);
                mu(i,j,k,l,:,:)  = cellextract(lens.mu_im);
                rho(i,j,k,l,:,:) = cellextract(lens.rho_im,'@(x)x(1,:)');
                % scale time delay
                y = abs(lens.loc_s(:,2));
                % Total time delay relative to optical axis
                tau = tau_im + repmat(0.5*y.^2,[1 size(tau_im,2)]);
                % Scaling
                t = tdelay(dist_mode,tau,lens.loc_s(:,2),zeta_0);
                % Plot (de)magnified signal against time
                siz = size(tau_im,2);
                signal = ones(size(tau_im));            % normalised amplitude per pulse
                period = 1e-6;                          % microsecond pulsar
                t_old = linspace(0, period*len, len)';  % pulse at each sourceloc
                t_new = repmat(t_old,[1 siz]) + t;      % replicate to nr. of images
                
            end % for l
        end % for k
    end % for j
end % for i
clear i j k l

% Extract time delays signals to compare between models
figure('visible','on'); hold on;
hSLines = scatter(t_old,signal(:,1),36,colour,'Marker','*','DisplayName','source');
set(gca,'Position',[.05 .05 .9 .9]);
for i=1:len
    hILines(i,:) = plot(t_new(i,:),mu_im(i,:).*signal(i,:),'o:','Color',colour(i,:));
end; clear i % for
hCGroup = hggroup; set(hILines,'Parent',hCGroup);
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
set(hCGroup,'DisplayName','images per source location');
xlabel('Time (s)'); ylabel('Flux relative to source'); legend('Location','Best');
mtit(char(title_str));

% Plots per model
title_str(1) = {['Source transiting ',lens.names.title]};
title_str(2) = {[lens.names.paramstr]};%,'N = ',num2str(num),', R = ',num2str(max_radius)]};
mtit(char(title_str));
end % function