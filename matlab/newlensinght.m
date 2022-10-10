% newlensing.m
%{
% function lens = newlensinght(sourceloc,dist_mode,params,num_model,varargin)
% flag_plot:
% 0: none,
% 1: 1 fig per source, legend, colourbar
% 2: 3 subplot per source: time delay surface and cross section,
magnification factor cross-section
% lens types given by num_model:
% 1: point mass lens, parameters: M
% 2: homogeneous disc, parameters: M, r_max / x_max
% 3: NFW lens, parameters: M, r_max / x_max, r_s / c
% flag_plot    = varargin{1};
    flag_summary = varargin{2};
    flag_pulsar  = varargin{3};
    scale_mode   = varargin{4};
%}

function lens = newlensinght(sourceloc,dist_mode,params,num_model,varargin)

global c G pc t_H eps_  % mks constants; single precision machine eps

if(nargin == 4)
    flag_plot    = 0;  % plot time delay surf, cross-sections
    flag_summary = 0;  % plot scaled tau(x), mu(x), rho(x) as source moves
    flag_pulsar  = 1;  % plot effect of lensing on pulsar signal
    scale_mode   = 1;  % if 1, use scaled lens params; if 0 use mks units
elseif(nargin >= 5)
    flag_plot    = varargin{1};
    flag_summary = varargin{2};
    flag_pulsar  = varargin{3};
    scale_mode   = varargin{4};
end % if

allcolors = get(0,'defaultAxesColorOrder');
allcolors(1,3) = 0.75; allcolors(3,1) = 0.75; % change 'red', 'blue'

% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
t_H = 3.09E17/0.6; % ! h=0.6 in Schneider
% machine epsilon
eps_ = double(eps('single'));

%% subplot setup
len = size(sourceloc,1); % nr of plots
if(flag_plot)
    if(flag_plot==1||len==2)
        nrows=1; ncols=3;  % 3 subplots
    elseif(flag_plot==2)
        nrows = ceil(len/4); ncols = 4;
    end
    ax = plot_axes(len,nrows,ncols);
end

%% Lens and model parameters

% Distances in m between source, lens, observer
switch(dist_mode{1})
    case('z')           % redshifts to source and lens
        Dd  = (c*t_H)*distang([0.27 0.73 0],0,dist_mode{2});
        Ds  = (c*t_H)*distang([0.27 0.73 0],0,dist_mode{3});
        Dds = (c*t_H)*distang([0.27 0.73 0],dist_mode{2},dist_mode{3});
        zd  = dist_mode{2};
    case('dist');       % distances in pc
        Dd  = dist_mode{2}*pc;
        Ds  = dist_mode{3}*pc;
        Dds = (dist_mode{3} - dist_mode{2})*pc;
        zd  = dist_mode{2}/(c*t_H); % assume d << D_H
end % switch
sigma_cr = c^2/(4*pi*G) * Ds/(Dd*Dds); % kg/m^2

if(isscalar(params)) % if only a mass is supplied, find other params
    
    %theta_0 = .009; % ang. radius of Einstein ring in arcsec.
    %M = theta_0^2 * (Dd*Ds)/Dds * c^2/(4*G);
    M = params(1)*M_earth; %3e-2*M_sol; % total lens mass
    r_E = sqrt(4*G) * sqrt(M/c^2) * sqrt(Dd*Dds/Ds);  % avoid catastrophic cancellation
    
    switch(num_model)
        case{1}
            zeta_0 = r_E;
        case{2}
            if(~scale_mode)
                r_max = 10^(-3)*pc;      % lens radius
                sigma_0 = M/(pi*(r_max)^2); % const. density in R3
                x_max = r_max/r_E
                kappa_0 = sigma_0 / sigma_cr;
            else
                x_max = 10             % r_max/r_E
                kappa_0 = x_max^(-2)    % scaled density
                % kappa_0 and x_max are NOT independent!
            end % if
            zeta_0 = r_E;
            
        case{3} % NFW
            %             % assume z << 1 so rho_crit equal to present-day value
            %             rho_crit = 3*t_H^(-2)/(8*pi*G); % ~ 9e-27 kg/m^3 (Hobson)
            %             % "virial radius" r_200 (nfw-cdm-halos p5)
            %             r_200 = (M/(200*rho_crit*4*pi/3))^(1/3);
            % r_s = 10^(-4.5)*pc;        % scale radius of lens (pc)
            %             if(~scale_mode)
            %                 % scale radius r_s, concentration c
            %                 x_max = 10^(-3)*pc/r_s % extent of lens (units of r_s)
            %                 c_nfw = r_200/r_s
            if(~scale_mode)
                sigma_cr = c^2./(4*pi*G) .* Ds./(Dd.*Dds); % kg/m^2
                r_0 = 10^(-4.5)*pc; % scale radius of lens (pc)
                r_max = 10^(-3)*pc; c_nfw = r_max/r_0; x_max = c_nfw;
                delta_c = M/((4*pi*r_0^3).*(log(1 + c_nfw) - c_nfw./(1+c_nfw)));
                kappa_s = delta_c*r_0/sigma_cr;
            end
            zeta_0 = r_0;
            clear delta_c rho_crit r_s
    end % switch
    
else % if params supplied, read them in
    M = params(1); zeta_0 = params(2);
    switch(num_model)
        case{1}
            % no other params
            
        case{2}
            x_max   = params(3);  % max. radius of lens
            kappa_0 = x_max^(-2); % scaled density
            % kappa_0 and x_max are NOT independent!
            
        case{3} % NFW
            x_max   = params(3); % max. radius of lens
            kappa_s = params(4); % Bartelmann density scaling param.
    end % switch
end % if

% Assign model details for plotting
M_str = num2str(M./M_sol,'%-0.3e');         % string for plotting
switch(num_model)
    case{1}
        lens.names = struct('scale','Einstein radii','title','Schwarzschild lens', ...
            'filestr','swzschild','paramstr',['Mass ',M_str,' M_{Sol}']);
        
    case{2}
        lens.names = struct('scale','Einstein radii','title','Disc profile', ...
            'filestr',strcat('hdisc_x_',stripdec(x_max,'.')), ...
            'paramstr',['Mass ',M_str,' M_{Sol}',' x_0: ',num2str(x_max*zeta_0/pc,'%-0.3e')]);
        
    case{3} % NFW
        lens.names = struct('scale','Scale radii r_{200}/c', ...
            'title','NFW profile', ...
            'filestr',strcat('nfw_x_',stripdec(x_max,'.'),'_kappa_',stripdec(kappa_s,',')), ...
            'paramstr',['M ',M_str,' M_{Sol} r_0 ',num2str(x_max*zeta_0/pc,'%-0.3e'),' pc \kappa_s ',num2str(kappa_s,'%-0.3e')]);
end % switch

% If scale_mode is off, scale source locations
if(~scale_mode); sourceloc = sourceloc./(zeta_0*Ds/Dd); end; % if
% Title for all plots incl. params and dist
title_str(1) = {['Source transiting ',lens.names.title]};
title_str(2) = {lens.names.paramstr};
title_str(3) = {['D_d = ',num2str(Dd/pc,'%-0.3e'),' pc, D_s = ',num2str(Ds/pc,'%-0.3e'),' pc']};

%% preallocate struct fields
[lens.loc_s(:,1), lens.loc_s(:,2)] = cart2pol(sourceloc(:,1),sourceloc(:,2));
lens.loc_im = cell(len,1); % numerical image loc
lens.tau_im = cell(len,1); % numerical tau at images
lens.mu_im  = cell(len,1); % numerical mu at images

switch(num_model)
    case{1}
        lens.rho_an = cell(len,1); % analytical image loc
        lens.tau_an = cell(len,1); % analytical tau at images
        lens.mu_an  = cell(len,1); % analytical mu at images
    case{2}
        lens.rho_an = cell(len,1); % analytical image loc
        lens.mu_an  = cell(len,1); % analytical mu at images
end % switch

%% Set up lens plane and convolution
if(exist('x_max','var')); % kappa(x) = 0 for x > x_max
    max_radius = max(x_max+1,2*ceil(max(abs(lens.loc_s(:,2)))));
else                      % kappa(x)=dirac(x), 100 arbitrary
    max_radius = max(100,2*ceil(max(abs(lens.loc_s(:,2)))));
end % if;
pix = 0.1; % max_radius / num % pixel size in scaled units (max. ~0.3)
num = ceil(max_radius/pix); % num pixels in [0,max_radius]
% Need to set R,N so that HT works for each SMD:
% kappa(R) = 0 and HT[kappa(R)] = sigma_fft(V) = 0
H = hankel_matrix(0,max_radius,num);
rho = H.r(:); x = [-flipud(H.r(:)); H.r(:)];
if(flag_plot); [x1,x2] = ndgrid(x,x); end; % if plotting plane, make grid

%% Numerical solution
% gravitational time delay
% 1. kernel is ln(rho): same dim as zero-padded potential
kernel = log(rho) - log(max(rho)); % lim(r-->inf) ker(r) --> 0
% 2. surface mass distribution: mass per pixel
kappa = kappa_fn(rho,num_model);
% 3. convolution
ht  = @(f) (H.C*(f(:)./H.m1)).*H.m2;
iht = @(F) (H.C*(F(:)./H.m2)).*H.m1;

kappa_fft = ht(kappa);
if(num_model==1) % kappa(r) = dirac(r) in cylindrical co-ords (projecting = \int dz)
    kappa_fft = 1/(2*pi)*ones(num,1); % manual HT; otherwise HT[r=0] NaN
end % if
krn_fft   = ht(kernel);
psi_fft   = 2*pi*kappa_fft .* krn_fft; % H[f1**f2] = 2pi.F1.F2
psi       = real ( iht( psi_fft ) ) ;
% adjust for krn(r) = ln(r) instead of ln(r-R)
psi       = psi + 2*pi*log(max(rho))*real(kappa_fft(1)); 

% rotate to get 2d potential time delay
if(flag_plot);
    [polyrev.X,polyrev.Y,polyrev.Z] = polyrevolve(rho,psi(1:num),0.4);
    polyrev.F = TriScatteredInterp(polyrev.X,polyrev.Y,polyrev.Z);
    td_pot    = polyrev.F(x1,x2);
end;

for i = 1:len
    phi_s = lens.loc_s(i,1); rho_s = lens.loc_s(i,2);
    sourcestr = sprintf('%0.5e %0.5e',sourceloc(i,:));
    
    % geometric time delay
    if(flag_plot); td_geom = 0.5.*((x1-sourceloc(i,1)).^2 + (x2-sourceloc(i,2)).^2); end;
    tau(:,3) = 0.5.*(x - rho_s).^2;
    
    % total time delay
    tau(:,2) = [flipud(psi(:));psi(:)];
    if(flag_plot); td_tot = td_geom - td_pot; end;
    tau(:,1) = tau(:,3) - tau(:,2); % can be < 0 because this is up to a const.
    % magnification factor
    mu = mu_fn(x);
    
    % image locations:
    % 1. mass distribution m(|x|) in subfunctions at end:
    
    % 2. Solve lens equation 
    x0 = rootsearch(@y_fn,@dydx_fn,@d2ydx2_fn,min(x),max(x))';
    
    % 3. From [x,phi_s] co-ords to polar co-ords [phi,rho]
    loc_im{i}(:,2) = abs(x0);
    loc_im{i}(:,1) = phi_s*ones(size(x0));
    ix = x0 < 0; % when x = -|x| = -r, rotate phi by pi
    loc_im{i}(ix,1) = phi_s - pi;  % set phi = phi_s - pi;
    clear ix

    % 4. interpolated solutions at image locations
    tau_im{i} = interp1(x,tau(:,1),x0,'cubic')';
    mu_im{i}  = interp1(x,mu,x0,'cubic')';
    
    %% Analytical soln.
    % radial co-ord rho; source pos. y
    y = rho_s;
    
    switch(num_model)
        case{1}
            phi_an = phi_s*ones(2,1); lens.rho_an{i} = zeros(2,1);
            x_an(1) = (y + sqrt( y^2 + 4 ) )/2 ;
            x_an(2) = (y - sqrt( y^2 + 4 ) )/2 ;
            lens.rho_an{i} = x_an;
            lens.tau_an{i} = 0.5*(x_an - y).^2 - log(abs(x_an));% + log(abs(max(x)));
            lens.mu_an{i}  = abs((1 - x_an.^(-4)).^(-1));
            
        case{2}
            if(x_max < 1)
                if(y < (1 - x_max^2)/x_max)
                    disp('y < (1 - x_max^2)/x_max')
                    phi_an  = phi_s*ones(3,1); lens.rho_an{i} = zeros(3,1);
                    x_an(1) = y * x_max^2/(x_max^2 - 1) ;
                    x_an(2) = (y + sqrt( y^2 + 4 ) )/2 ;
                    x_an(3) = (y - sqrt( y^2 + 4 ) )/2 ;
                elseif(y == (1 - x_max^2)/x_max)
                    phi_an  = phi_s*ones(2,1); lens.rho_an{i} = zeros(2,1);
                    x_an(1) = x_max ;
                    x_an(2) = 1/x_max ;
                elseif(y > (1 - x_max^2)/x_max)
                    phi_an  = phi_s;
                    x_an    = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            elseif(x_max > 1)
                phi_an = phi_s;
                if(y <= (x_max^2 - 1)/x_max)
                    x_an = y * x_max^2/(x_max^2 - 1) ;
                else
                    x_an = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            else % if x_max = 1 all points are focused onto the optical axis
                phi_an = phi_s; x_an = 0;
            end
            % SMD kappa = surface mass per pix = surface density * pix. area
            % So x_0^(-2) = kappa_0 --> kappa_0*pix^2
            mu_an = ((1 - 1./x_an.^2).*(1 + 1./x_an.^2)).^(-1); % true x_im > x0
            [lens.rho_an{i},ix] = sort(x_an);
            lens.mu_an{i} = abs(mu_an(ix));
            indx = abs(lens.rho_an{i}) < x_max;  % logical indexing
            lens.mu_an{i}(indx) = (1 - kappa_0).^(-2); % for x_im < x0
    end
    
    %% Plot results
    if(flag_plot)
        %% time delay surface
        %%{
        if(flag_plot==1)
            if(i==1)
                figure('visible','on');
                title(strcat('Time delay surface',lens.names.title));
                colorbar('location','southoutside');
                legend('Location','NorthOutside');
                xlabel(lens.names.scale); ylabel(lens.names.scale);
            end; % if
            subplot(nrows,ncols,1);
        elseif(flag_plot==2)
            figure(1);
            subplot(nrows,ncols,i);
            set(gca,'pos',ax.coord(i,:));
        end % if
        hold on;
        hTimeDelay = pcolor(x1,x2,td_tot);
        colormap(flipud(thermal)); shading interp;
        set(hTimeDelay,'DisplayName','Time delay surface');
        
        % plot numerical (x) image locations
        [loc_x,loc_y] = pol2cart(loc_im{i}(:,1),loc_im{i}(:,2));
        hSLines = scatter(loc_x,loc_y,[],allcolors(1:length(loc_x),:), ...
            'Marker','x','DisplayName','Numerical image location');
        clear loc_x loc_y
        
        % plot analytical (o) image locations
        if(isfield(lens,'rho_an'))
            [loc_x,loc_y] = pol2cart(phi_an(:),lens.rho_an{i}(:));
            hSLines = scatter(loc_x,loc_y,[],allcolors(1:length(loc_x),:), ...
                'Marker','o','DisplayName','Analytical image location');
            clear loc_x loc_y
        end
        % if the lens has finite radius, plot lens
        if(exist('x_max','var')) % if the lens has finite radius
            hDisc = polar(linspace(0,2*pi,500)',x_max*ones(500,1),'k-');
            set(hDisc,'DisplayName','extent of lens');
        end
        % plot source
        hSource = polar(phi_s,rho_s,'ks');
        switch(flag_plot)
            case{1}
                set(hSource,'DisplayName',strcat('source at',sourcestr));
            case{2}
                set(hSource,'DisplayName','source');
        end;
        clear h*; hold off; axis image;
        if(i==len); colorbar('Location','NorthOutside'); end; % if
        %}
        
        %% time delay cross-section
        switch(flag_plot)
            case{1}
                subplot(nrows,ncols,2);
                if(i==1)
                    title('Radial cross-section of time delay surface');
                    xlabel(lens.names.scale); ylabel('Scaled time delay');
                end; % if
            case{2}
                figure(2); subplot(nrows,ncols,i);
                set(gca,'pos',ax.coord(i,:));
        end % switch
        hold on;
        % time delay cross-section
        plot(x,tau(:,1),'k-', x,tau(:,2),'k-.', x,tau(:,3),'k:');
        for l=1:size(tau_im{i},1)
            hSLines(l) = plot([-num/2*pix num/2*pix], ...
                [tau_im{i}(l) tau_im{i}(l)], '--','color',allcolors(l,:));
        end; clear l
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        hold off; axis([-15 15 0 200]); clear h*
        if(flag_plot==1)
            legend('Total time delay','Fermat potential','Geometric contribution', ...
                'Delay at images','Location','Best');
        end; % if
        
        %% magnification factor cross-section
        switch(flag_plot)
            case{1}
                subplot(nrows,ncols,3);
                if(i==1)
                    title('Radial cross-section of magnification surface');
                    xlabel(lens.names.scale); ylabel('Scaled magnification factor');
                    legend('Location','BestOutside');
                end; % if
            case{2}
                figure(3); subplot(nrows,ncols,i);
                set(gca,'pos',ax.coord(i,:));
        end % switch
        
        hold on;
        plot(x,mu,'k-','DisplayName','Magnification factor');
        
        % plot numerical (x) image magnification
        for l=1:size(loc_im{i},2)
            hSLines(l) = plot(loc_im{i}(2,l),mu_im{i}(l),'x','color',allcolors(l,:));
        end; clear l % for
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(hSGroup,'DisplayName','Numerical magnification factor');
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        
        % plot analytical (o) image locations
        if(isfield(lens,'rho_an'))
            for l=1:size(lens.rho_an{i},1)
                hCLines(l) = plot(lens.rho_an{i}(l),lens.mu_an{i}(l),'o','color',allcolors(l,:));
            end; clear l % for
            hCGroup = hggroup; set(hCLines,'Parent',hCGroup);
            set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(hCGroup,'DisplayName','Analytical image location');
        end % if
        hold off; axis tight; clear h*
        
        %         % Save plot
        %         if(flag_plot==1) % new plot for each source
        %             source_cell = stripdec(source,'.');
        %             [A B] = deal(source_cell{:});
        %             sstr = cat(2,'(',A,' ',B,')');
        %             filename = strcat('plot_',lens.names.filestr,'_',sstr);
        %             saveas(gcf,filename,'png');
        %         end % if
        clear mu tau
    end % if
    
end; clear i % for

%% Plot brightness and time delay curves during transit

% plot does not work with cell types
lens.tau_im = cellextract(tau_im);
lens.mu_im  = cellextract(mu_im);
lens.x_im   = cellextract(loc_im,'@(x)x(:,2)');
lens.phi_im = cellextract(loc_im,'@(x)x(:,1)');
clear *_im
% plotting
if(flag_summary)
    plot_summary(lens,dist_mode);
    mtit(char(title_str)); % for A \n B C use char(A,[B,C]);
    clear title_str
end % if

%% Impose time delay on signal
if(flag_pulsar)
    % Total time delay relative to optical axis
    clear tau;
    tau = lens.tau_im + repmat(0.5*abs(lens.loc_s(:,2)).^2, [1 size(lens.tau_im,2)]);
    % Scaling
    t = zeta_0^2/c * Ds/(Dd*Dds) * (1+zd).*tau;
    % Plotting
    plot_tdelay({t},{lens.mu_im},hypot(sourceloc(:,1),sourceloc(:,2)));
    dist_str = mat2str(cell2mat(dist_mode(2:end-1)),3)/1e3;
    title_str(1) = {[lens.names.title,lens.names.paramstr]};
    title_str(2) = {['Lens dist: ',dist_str,' kpc; Source dist: ',num2str(Ds/(1e3*pc)),' kpc']};
    title_str = regexprep(title_str,'e\+0',' \\times 10^');
    title(gca,title_str);
    
end % if

%% Surface mass density
    function kappa = kappa_fn(rho,num_model)
        kappa = zeros(size(rho)); % convergence (SMD)
        
        switch(num_model)
            case{1}
                kappa(1) = pi; kappa(num) = 0;
                
            case{2}
                % kappa(x) = sigma(zeta_0 * x)/sigma_cr NOT sigma(x)/sigma_cr
                kappa(rho <= x_max) = kappa_0;% * pix^2; % mass = den x area
                % indexing returns row vect.: reshape
                kappa = reshape(kappa,size(rho));
                
            case{3}
                % conditional eqn. for f (Bartelmann)
                F = zeros(size(rho)); % surface density is 0, x > x_max
                ind = find((rho > 1) & (rho <= x_max));
                F(ind) = 1 - 2./sqrt(rho(ind).^2 - 1).*atan(sqrt((rho(ind) - 1)./(1 + rho(ind))));
                F(rho < 1) = 1 - 2./sqrt(1 - rho(rho < 1).^2).*atanh(sqrt((1 - rho(rho < 1))./(1 + rho(rho < 1))));
                F(rho ==1) = 0;
                kappa = 2*kappa_s.*F./(rho.^2 - 1);
                
        end % switch
    end % function

%% Magnification factor
    function mu_x = mu_fn(x)
        X = abs(x); m = ones(size(x)); kappa_x = zeros(size(x));
        
        switch(num_model)
            case{1}
                if(~isempty(X(X < eps_)))
                    kappa_x(X < eps_) = 1; % zero except at x=0
                end; % if
                
            case{2} % m/x = x/x_max^2; 1/x
                m(X <= x_max)       = x(X < x_max).^2./x_max.^2;
                kappa_x(X <= x_max) = x_max^(-2); % kappa_0 * pix^2;
                
            case{3}
                ind = find((X > 1) & (X <= x_max));
                g = zeros(size(X));
                g(ind)   = 2./sqrt(X(ind).^2 - 1) .* atan(sqrt((X(ind) - 1)./(X(ind) + 1)));
                g(X < 1) = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                g(X ==1) = 1;
                m(X <= x_max) = 4*kappa_s*(log(X(X <= x_max)/2) + g(X <= x_max)); 
                m(X==0)=0; 

                f = zeros(size(X));
                f(ind) = 1 - 2./sqrt(X(ind).^2 - 1).*atan(sqrt((X(ind) - 1)./(1 + X(ind))));
                f(X < 1) = 1 - 2./sqrt(1 - X(X < 1).^2).*atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                f(X ==1) = 0;
                kappa_x(X <= x_max)  = 2*kappa_s.*f(X <= x_max)./(X(X <= x_max).^2 - 1);
        end
        detA = (1 - m./x.^2).*(1 + m./x.^2 - 2*kappa_x); % eqn. 8.16 SE&F
        mu_x = abs(detA.^(-1)); % |inverse| = magnification
    end % function

%% Functions to define y, dydx, d2ydx2

    function y = y_fn(x)
        X = abs(x);
        m = ones(size(x));
        switch(num_model)
            case{1}
                % no changes
                
            case{2} % m/x = x/x_max^2; 1/x
                m(X <= x_max) = X(X <= x_max).^2./x_max.^2;
                
            case{3}
                ind = find((X > 1) & (X <= x_max));
                g = zeros(size(X));
                g(ind)   = 2./sqrt(X(ind).^2 - 1) .* atan(sqrt((X(ind) - 1)./(X(ind) + 1)));
                g(X < 1) = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                g(X ==1) = 1;
                m(X <= x_max) = 4*kappa_s*(log(X(X <= x_max)/2) + g(X <= x_max)); 
                m(X==0)=0; 
        end
        % rho_s = f(x) same as 0 = f(x) - rho_s
        y = x - m./x - rho_s; 
        if(num_model==3); y(x==0) = -rho_s; end; % limit
    end % function

    function dydx = dydx_fn(x)
        X = abs(x);
        m = ones(size(x)); dmdx = zeros(size(x));
        switch(num_model)
            case{1}
                % no changes necessary
                
            case{2}
                m(X <= x_max) = X(X <= x_max).^2./x_max.^2;
                dmdx(X <= x_max) = 2.*X(X <= x_max)./x_max^2;
                
            case{3}
                g = zeros(size(X)); dg=g; 
                ind      = find((X > 1) & (X <= x_max));
                g(ind)   = 2./sqrt(X(ind).^2-1) .* atan(sqrt((X(ind)-1)./(X(ind)+1)));
                dg(ind)  = -2.*(sqrt((X(ind)-1)./(X(ind)+1)) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1))) + X(ind)-1) ./ ((X(ind)-1).^2.*X(ind).*(X(ind)+1));
                clear ind
                ind      = find(X < 1);
                g(ind)   = 2./sqrt(1-X(ind).^2) .* atanh(sqrt((1-X(ind))./(1+X(ind))));
                dg(ind)  = 2.*(sqrt((1-X(ind))./(X(ind)+1)) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) + X(ind)-1) ./ (X(ind).*sqrt((1-X(ind))./(X(ind)+1)) .* (1-X(ind).^2).^(3/2));
                clear ind
                g(X==1)  = 1;
                dg(X==1) = 0;
                
                m(X <= x_max)    = 4*kappa_s*(log(X(X <= x_max)./2) + g(X <= x_max)); m(X==0)=0;
                dmdx(X <= x_max) = 4*kappa_s*(1./(X(X <= x_max)./2) + dg(X <= x_max)); dmdx(X==0)=0;
        end
        dydx = 1 - (x.*dmdx - m)./x.^2; % diverges at x=0
    end % function

    function d2ydx2 = d2ydx2_fn(x)
        X = abs(x); m = ones(size(x)); 
        dmdx = zeros(size(x)); d2mdx2 = zeros(size(x));
        switch(num_model)
            case{1}
                % no changes
                
            case{2}
                m(X <= x_max)      = X(X <= x_max).^2./x_max.^2;
                dmdx(X <= x_max)   = 2.*X(X <= x_max)./x_max^2;
                d2mdx2(X <= x_max) = 2/x_max^2;
                
            case{3}
                g = zeros(size(X)); dg=g; d2g=g;
                ind       = find((X > 1) & (X <= x_max));
                g(ind)    = 2./sqrt(X(ind).^2-1) .* atan(sqrt((X(ind)-1)./(X(ind)+1)));
                dg(ind)   = -2*(sqrt((X(ind)-1)./(X(ind)+1)) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1)))+X(ind)-1) ./ ((X(ind)-1).^2.*X(ind).*(X(ind)+1));
                d2g(ind)  = (-4*X(ind).^3 + 4*X(ind).^2 + 2*sqrt((X(ind)-1)./(X(ind)+1)) .* (2*X(ind).^2+1) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1))) +X(ind)-1) ./ (X(ind).^2 .* sqrt((X(ind)-1)./(X(ind)+1)) .* (X(ind).^2-1).^(5/2));
                clear ind
                ind       = find(X < 1);
                g(ind)    = 2./sqrt(1-X(ind).^2) .* atanh(sqrt((1-X(ind))./(1+X(ind))));
                dg(ind)   = 2.*(sqrt((1-X(ind))./(X(ind)+1)) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) + X(ind)-1) ./ (X(ind).*sqrt((1-X(ind))./(X(ind) + 1)) .* (1-X(ind).^2).^(3/2));
                d2g(ind)  = (4.*X(ind).^3 - 4.*X(ind).^2 + 2.*sqrt((1-X(ind))./(X(ind)+1)) .* (2*X(ind).^2+1) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) - X(ind)+1) ./ (X(ind).^2 .* sqrt((1-X(ind))./(X(ind)+1)) .* (1-X(ind).^2).^(5/2));
                clear ind
                g(X == 1) = 1; dg(X == 1) = 0; d2g(X == 1) = 0;
                
                m(X <= x_max) = 4*kappa_s*(log(X(X <= x_max)./2) + g(X <= x_max)); 
                dmdx(X <= x_max) = 4*kappa_s*(2./X(X <= x_max) + dg(X <= x_max)); 
                d2mdx2(X <= x_max) = 4*kappa_s*(-2./X(X <= x_max).^2 + d2g(X <= x_max)); 
                m(X==0)=0; dmdx(X==0)=0; % d2mdx2 diverges at x=0
        end
        d2ydx2 = (-2.*m + 2.*x.*dmdx - d2mdx2)./x.^3; % diverges at x=0
        if(num_model==3); d2ydx2(X==0) = -Inf; end;
    end % function

end % function