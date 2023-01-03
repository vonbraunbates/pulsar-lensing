% newlensing.m
%{
% tim = newlensing(sourceloc,num_model,flag_plot)
% flag_plot:
% 0: none,
% 1: 1 fig per source, legend, colourbar
% 2: 3 subplot per source: time delay surface and cross section,
magnification factor cross-section
% lens types given by num_model:
% 1: point mass lens, parameters: M
% 2: homogeneous disc, parameters: M, r_max / x_max
% 3: NFW lens, parameters: M, r_max / x_max, r_s / c
% scale_mode: optional arg
% .true.: use scaled quantities (default)
% .false.: use physical quantites
%}

function lens = newlensing(sourceloc,num_model,flag_plot,scale_mode)

global c G pc t_H eps_ % mks universal constants; single precision machine eps

if(nargin==3) scale_mode=1; end;
allcolors = get(0,'defaultAxesColorOrder');
allcolors(1,3) = 0.75; allcolors(3,1) = 0.75; % change 'red', 'blue'

% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
t_H = 3.09E17/0.6; % ! h=0.6 in Schneider
% machine epsilon 
eps_ = double(eps('single'));

%% subplot setup
if(flag_plot)
    len = size(sourceloc,1); % nr of plots
    if(flag_plot==1||len==2)
        nrows=1; ncols=3;  % 3 subplots
    elseif(flag_plot==2)
        nrows = ceil(len/6); ncols = 6;
    end
    % set figure size, visibility
    set(0,'DefaultFigureVisible','off');
    scrsz = get(0,'ScreenSize');
    set(0,'DefaultFigurePosition',[1 .9*scrsz(4) .99*scrsz(3) .9*scrsz(4)]);
    % axes in plot
    ax.min = 0.05; ax.max = 0.95; ax.gap = 0.05;
    ax.size = (ax.max - ax.min)./[ncols nrows];
    ax.box = ax.size - ax.gap;
    ax.coord(:,1) = ax.min + ax.size(1).*mod([1:len]-1,ncols); % x starting co-ord
    ax.coord(:,2) = ax.max - ax.size(2).*ceil([1:len]./ncols); % y starting co-ord
    ax.coord(:,3) = ax.box(1); % width
    ax.coord(:,4) = ax.box(2); % height
    ax.coord(end+1,:) = [0 0 1 1]; % figure axes

    % Remove warnings in legend
    warning('off','MATLAB:legend:UnsupportedFaceColor');
    warning('off','MATLAB:legend:PlotEmpty');
    warning('off','MATLAB:legend:IgnoringExtraEntries');
end

%% Lens plane and model setup
% distances to source and lens
Ds = 2e3*pc; Dd = 1e3*pc; % m
Dds = Ds - Dd; % assume d << D_H
sigma_cr = c^2/(4*pi*G) * Ds/(Dd*Dds); % kg/m^2

% lens, source parameters
M = 1e4*M_earth; % total lens mass
mstr = num2str(M/M_sol,'%6.4e'); % string for plotting
r_E = sqrt(4*G*M/c^2*Dd*Dds/Ds); % Einstein radius
zeta_0 = r_E; % lens plane: x = zeta/zeta_0

switch(num_model)
    case{1}
        model = struct('scale','Einstein radii','title','Schwarzschild lens', ...
            'filestr','swzschild','paramstr',['Mass ',mstr,' M_{\odot}']);
        params = struct('mass',M,'r_E',r_E);
        
    case{2}
        if(~scale_mode)
            r_max = r_E; 
            sigma_0 = M/(pi*(r_max)^2); % const. density
            x_max = r_max/zeta_0;
            kappa_0 = sigma_0 / sigma_cr;
        else
            x_max = 0.9 % r_max/zeta_0
            kappa_0 = x_max^(-2) % scaled density
            % kappa_0 and x_max are NOT independent!
        end % if
        model = struct('scale','Einstein radii','title','Homogeneous disc', ...
            'filestr',strcat('hdisc_x_',stripdec(x_max,'.')), ...
            'paramstr',['Mass ',mstr,' M_{\odot}',' x_0: ',num2str(x_max,'%6.4e')]);
        params = struct('mass',M,'r_E',r_E,'x_max',x_max);
        
    case{3} % NFW
        % assume z << 1 so rho_crit equal to present-day value
        rho_crit = 3*t_H^(-2)/(8*pi*G); % ~ 9e-27 kg/m^3 (Hobson)
        % "virial radius" r_200 (nfw-cdm-halos p5)
        r_200 = (M/(200*rho_crit*4*pi/3))^(1/3);
        if(~scale_mode)
            % scale radius r_s, concentration c
            r_s = 10^(-4.5)*pc; % cutoff radiusof lens
            x_max = r_s/zeta_0; % scaled radius of lens
            c_nfw = r_200/r_s;
        else % param. is c
            c_nfw = 2e6; % c > 1 (wright-nfw-halos p3)
            x_max = .5; % turnover radius r_max/r_s
        end
        % characteristic overdensity delta_c (nfw-cdm-halos p5)
        delta_c = (200/3) * c_nfw^3 / (log(1+c_nfw) - c_nfw/(1+c_nfw));
        kappa_s = delta_c*rho_crit*r_200/c/sigma_cr;
        model = struct('scale','Scale radii r_{200}/c', ...
            'title','Navarro-Frenck-White profile', ...
            'filestr',strcat('nfw_x_',stripdec(x_max,'.'),'_c_',stripdec(c_nfw,',')), ...
            'paramstr',['Mass ',mstr,' M_{\odot}',' turnover radius: ',num2str(x_max,'%6.4g'),' c: ',num2str(c_nfw,'%6.4g')]);
        params = struct('mass',M,'r_E',r_E,'x_max',x_max,'c',c);
        clear delta_c rho_crit r_s
end % switch

lens.params = params; % save model parameters
lens.names = model; % save model details
clear model params

% set up grid for lens plane
num = pow2(10); % unpadded box will be num x num px.
boxlim = 16; % padded box edge is [-xlim xlim]
pix = boxlim / num  % pixel size in scaled units
% zero-padded grid (for grav. time delay)
[x1_big,x2_big] = ndgrid(linspace(-boxlim,boxlim,2*num));
[~,rho_big] = cart2pol(x1_big,x2_big); % incl. padding (for kernel)
clear xlim
% remove padding (for geom. time delay)
x1 = x1_big(num/2+1:3*num/2,num/2+1:3*num/2);
x2 = x2_big(num/2+1:3*num/2,num/2+1:3*num/2);
rho = rho_big(num/2+1:3*num/2,num/2+1:3*num/2);
clear x*_big

% preallocate struct fields
[lens.loc_s(:,1), lens.loc_s(:,2)] = cart2pol(sourceloc(:,1),sourceloc(:,2));
lens.loc_im = cell(size(sourceloc,1),1);
lens.tau_im = cell(size(sourceloc,1),1);
lens.mu_im  = cell(size(sourceloc,1),1);
lens.tau    = cell(size(sourceloc,1),1);
lens.mu     = cell(size(sourceloc,1),1);
lens.rho_an = cell(size(sourceloc,1),1);
lens.tau_an = cell(size(sourceloc,1),1);
lens.mu_an  = cell(size(sourceloc,1),1);

%% Numerical solution

% gravitational time delay
% 1. kernel is ln(rho): same dim as zero-padded potential
kernel = log(rho_big);
kernel = circshift( kernel , [-num+1 -num+1] );

% 2. surface mass distribution: mass per pixel
kappa = zeros(size(rho)); % convergence (SMD) without zp (at this stage!)

switch(num_model)
    case{1}
        % \int_R dx dirac(x) = 1 iff dirac(0) = 1/(\int_R dx)
        % normally R is the central pixel, area 1/(pix)^2 BUT
        % we want independence from grid fineness so set area(R)=1
        dirac_x = 1;
        % dirac(at) = |a|^(-n) * dirac(t)
        dirac_zeta = zeta_0^(-2) * dirac_x;
        % kappa = sigma(zeta)/sigma_cr = M*dirac(zeta)/sigma_cr
        % dirac(0) in centre of rho NOT rho_big
        kappa(num/2,num/2) = M*dirac_zeta/sigma_cr;
        
    case{2}
        % kappa(x) = sigma(zeta_0 * x)/sigma_cr NOT sigma(x)/sigma_cr
        kappa(rho <= x_max) = kappa_0 * pix^2; % mass = den x area
        % indexing returns row vect.: reshape
        kappa = reshape(kappa,size(rho));
        
    case{3}
        % conditional eqn. for f (Bartelmann)
        F(rho > 1) = 1 - 2./sqrt(rho(rho > 1).^2 - 1).*atan(sqrt((rho(rho > 1) - 1)./(1 + rho(rho > 1))));
        F(rho < 1) = 1 - 2./sqrt(1 - rho(rho < 1).^2).*atanh(sqrt((1 - rho(rho < 1))./(1 + rho(rho < 1))));
        F(rho ==1) = 0;
        % indexing returns row vect.: reshape
        F = reshape(F,size(rho)); 
        kappa = 2*kappa_s.*F./(rho.^2 - 1);
        
end % switch

kappa(2*num,2*num) = 0; % zero-padding

% 3. integrate SMD to get potential: equal to convolution in k-space
% F(f)F(g) = F(f*g) = F(\int dxdy f(x,y)g(x,y)) then invert to get \int
sigma_fft = fft2( kappa ) ; krn_fft = fft2( kernel ) ;
psi_fft = sigma_fft .* krn_fft ; % F[f**g] = 2pi*F[f]F[g]
psi = real ( ifft2( psi_fft /pi) ) ; % inverse transform multiplies by pi
clear kernel *fft
% remove padding
td_pot = psi(1:num,1:num); 
clear psi*

% 1. minimise calculations by using smallest possible subset of R2
rho_data = rho(logical(eye(num)));
rho_data = [(-rho_data(1:num/2));rho_data(num/2+1:num)];
x = rho_data;

for i = 1:size(sourceloc,1)
    source = sourceloc(i,:); % dimensionless source loc.
    sourcestr = sprintf('%0.5e %0.5e',source);
    [phi_s rho_s] = cart2pol(source(1),source(2)); % source in polar co-ords
    if(phi_s < 0); lens.loc_s(i,2) = -rho_s; end; 
    
    % geometric time delay
    td_geom = 0.5.*((x1-source(1)).^2 + (x2-source(2)).^2);
    
    % total time delay
    td_tot = td_geom - td_pot;
    
    % image locations:
    [xi,yi] = pol2cart(phi_s*ones(size(rho_data)),rho_data);
    % 1. mass distribution m(|x|)
    %{
    This goes in subfunctions at end:
    y = @(x)(x - m/x);
    dydx = @(x)(1 - (x*dmdx - m)/x^2);
    d2ydx2 = @(x)((-2*m + 2*x*dmdx - d2mdx2)/x^3);
    mu = @(x)([(1 - m/x^2)(1 + m/x^2 - 2*kappa)]^(-1));
    
    % 3. Solve lens equation NUMERICALLY
    generic GL catastrophes lead to:
    cusps: 2 closely-spaced roots of f --> 1 root of f'
    folds: 3 closely-spaced roots of f --> 1 root of f"
    find roots of each derivative
    %}
    x0 = rootsearch(@y_fn,@dydx_fn,@d2ydx2_fn,min(x),max(x))';
    %x0 = reshape(x0,numel(x0),1); % col. vect
    
    % 3. reconvert to Cartesian co-ords
    lens.loc_im{i}(1,:) = x0; lens.loc_im{i}(2,:) = phi_s*ones(size(x0));
    [x_loc y_loc] = pol2cart(lens.loc_im{i}(2,:),lens.loc_im{i}(1,:));
    
    % 4. Cross-section at phi=phi_s for plotting
    lens.tau{i}(:,1) = interpn(x1,x2,td_tot,xi,yi);
    lens.tau{i}(:,2) = interpn(x1,x2,td_pot,xi,yi);
    lens.tau{i}(:,3) = interpn(x1,x2,td_geom,xi,yi);
    lens.mu{i} = mu_fn(x);
    
    % 5. interpolated solutions at image locations
    lens.tau_im{i} = interpn(x1,x2,td_tot,x_loc,y_loc);
    lens.mu_im{i} = interp1(x,lens.mu{i},x0,'spline');
    
%%{    
    %% Analytical soln.
    % radial co-ord rho; source pos. y
    y = rho_s;
    
    % locations of extrema
    switch(num_model)
        case{1}
            phi_im = phi_s*ones(2,1); rho_im = zeros(2,1);
            rho_im(1) = (y + sqrt( y^2 + 4 ) )/2 ;
            rho_im(2) = (y - sqrt( y^2 + 4 ) )/2 ;
            rho_im = sort(rho_im);
            [x_im,y_im] = pol2cart(phi_im,rho_im);
            lens.rho_an{i} = rho_im;
            lens.tau_an{i} = 0.5*(rho_im - y).^2 - log(abs(rho_im));
            lens.mu_an{i}  = abs((1 - rho_im.^(-4)).^(-1));
            clear *_im
            
        case{2}
            if(x_max < 1)
                if(y < (1 - x_max^2)/x_max)
                    disp('y < (1 - x_max^2)/x_max')
                    phi_im = phi_s*ones(3,1); rho_im = zeros(3,1);
                    rho_im(1) = y * x_max^2/(x_max^2 - 1) ;
                    rho_im(2) = (y + sqrt( y^2 + 4 ) )/2 ;
                    rho_im(3) = (y - sqrt( y^2 + 4 ) )/2 ;
                elseif(y == (1 - x_max^2)/x_max)
                    phi_im = phi_s*ones(2,1); rho_im = zeros(2,1);
                    rho_im(1) = x_max ;
                    rho_im(2) = 1/x_max ;
                elseif(y > (1 - x_max^2)/x_max)
                    phi_im = phi_s;
                    rho_im = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            elseif(x_max > 1)
                phi_im = phi_s;
                if(y <= (x_max^2 - 1)/x_max)
                    rho_im = y * x_max^2/(x_max^2 - 1) ;
                else
                    rho_im = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            else % if x_max = 1 all points are focused onto the optical axis
                phi_im = phi_s; rho_im = 0;
            end
            lens.rho_an{i} = sort(rho_im);
            [x_im,y_im] = pol2cart(phi_im,rho_im);
            lens.mu_an{i} = abs((1 - kappa_0).^(-2).*ones(size(rho_im)));
    end
    
    %% Plot results
    if(flag_plot)
        switch(flag_plot)
            case{1}
                figure(1); subplot(nrows,ncols,1);
                title(strcat('Time delay surface',lens.names.title));
                colorbar('location','southoutside');
                legend('Location','NorthOutside');
                xlabel(lens.names.scale); ylabel(lens.names.scale);
            case{2}
                figure(1);
                subplot(nrows,ncols,i);
                set(gca,'pos',ax.coord(i,:));
        end % switch
        
        hold on;
        %% time delay surface
        hTimeDelay = pcolor(x1,x2,td_tot);
        colormap(flipud(thermal)); shading interp;
        set(hTimeDelay,'DisplayName','Time delay surface');
        % plot numerical (x) image locations
        for l=1:size(lens.loc_im{i},1)
            hSLines(l) = plot(x_loc(l),y_loc(l),'x','color',allcolors(l,:));
        end; clear l
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(hSGroup,'DisplayName','Numerical image location');
        % plot analytical (o) image locations
        if(exist('rho_im','var'))
            for l=1:length(rho_im)
                hCLines(l) = plot(x_im(l),y_im(l),'o','color',allcolors(l,:));
            end; clear l
            hCGroup = hggroup; set(hCLines,'Parent',hCGroup);
            set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(hCGroup,'DisplayName','Analytical image location');
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
        clear h*; hold off; axis image; % tight axes
        
        %% time delay cross-section
        switch(flag_plot)
            case{1}
                figure(1); subplot(nrows,ncols,2);
                title('Radial cross-section of time delay surface');
                xlabel(lens.names.scale); ylabel('Scaled time delay');
            case{2}
                figure(2); subplot(nrows,ncols,i);
                set(gca,'pos',ax.coord(i,:));
        end % switch
        hold on;
        % time delay cross-section
        plot(rho_data,lens.tau{i}(:,1),'k-', ...
            rho_data,lens.tau{i}(:,2),'k-.', ...
            rho_data,lens.tau{i}(:,3),'k:');
        for l=1:size(lens.tau_im{i},1)
            hSLines(l) = plot([-num/2*pix num/2*pix], ...
                [lens.tau_im{i}(l) lens.tau_im{i}(l)], '--','color',allcolors(l,:));
        end; clear l
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        hold off; axis tight; clear h*
        if(flag_plot==1)
            legend('Total time delay','Fermat potential','Geometric contribution', ...
             'Delay at images','Location','Best');
        end; % if
        
        %% magnification factor cross-section
        switch(flag_plot)
            case{1}
                figure(1); subplot(nrows,ncols,3);
                title('Radial cross-section of magnification surface');
                xlabel(lens.names.scale); ylabel('Scaled magnification factor');
                legend('Location','BestOutside');
            case{2}
                figure(3); subplot(nrows,ncols,i);
                set(gca,'pos',ax.coord(i,:));
        end % switch
        
        hold on;
        plot(rho_data,lens.mu{i},'k-','DisplayName','Magnification factor');
        
        % plot numerical (x) image magnification
        for l=1:size(lens.loc_im{i},1)
            hSLines(l) = plot(lens.loc_im{i}(l,1),lens.mu_im{i}(l),'x','color',allcolors(l,:));
            if(exist('rho_im','var'))
                fprintf('analytical value: %.5f ',mu_im(l));
                fprintf('interpn 1d mu: %.5f ',lens.mu_im{i}(l));
                fprintf(['rel. difference: %.5f ','spline','\n'],100*abs(lens.mu_im{i}(l) - mu_im(l))/mu_im(l));
            end
        end; clear l % for
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(hSGroup,'DisplayName','Numerical magnification factor');
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        
        % plot analytical (o) image locations
        if(exist('rho_im','var'))
            for l=1:length(rho_im)
                hCLines(l) = plot(rho_im(l),mu_im(l,1),'o','color',allcolors(l,:));
                if(size(mu_im,2)==2)
                    if(abs(mu_im(l,2) - mu_im(l,1)) > eps_); disp('Analytical solution error!'); end;
                end % if
            end; clear l % for
            hCGroup = hggroup; set(hCLines,'Parent',hCGroup);
            set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(hCGroup,'DisplayName','Analytical image location');
        end % if
        hold off; axis tight; clear h*
%         
%         % Save plot
%         if(flag_plot==1) % new plot for each source
%             source_cell = stripdec(source,'.');
%             [A B] = deal(source_cell{:});
%             sstr = cat(2,'(',A,' ',B,')');
%             filename = strcat('plot_',lens.names.filestr,'_',sstr);
%             saveas(gcf,filename,'png');
%         end % if
%     end % if
%}
    end
end; clear i % for

%% Plot brightness and time delay curves during transit

% plot does not work with cell types
tau_im = cellextract(lens.tau_im);
mu_im = cellextract(lens.mu_im);
rho_im = cellextract(lens.loc_im,'@(x)x(1,:)');
phi_im = cellextract(lens.loc_im,'@(x)x(2,:)');

figure('visible','on');

% analytical plots, if soln extant
if(isfield(lens,'rho_an'))
    rho_an = cellextract(lens.rho_an);
    subplot(2,2,1); hold all;
    fLines = plot(lens.loc_s(:,2),rho_an,'-');
    %fLines = plot(lens.loc_s(:,2),(rho_an - rho_im)./rho_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('r(x_s)'), title('Fractional error in radial image locations');
end

if(isfield(lens,'tau_an'))
    tau_an = cellextract(lens.tau_an);
    subplot(2,2,2); hold all;
    %plot(lens.loc_s(:,2),0.5*(rho_an - repmat(lens.loc_s(:,2),1,2)).^2,'k--');
    fLines = plot(lens.loc_s(:,2),tau_an,'-');
    %fLines = plot(lens.loc_s(:,2),(tau_an - tau_im)./tau_an,'o');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Fractional error in time delay');
    %%{
    subplot(2,2,4); hold all;
    fLines = plot(lens.loc_s(:,2),(tau_an - tau_im)./tau_an,'o');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Absolute time delay');
    legend('Location','best'), axis tight;
    %}
end

if(isfield(lens,'mu_an'))
    mu_an = cellextract(lens.mu_an);
    subplot(2,2,3); hold all;
    fLines = plot(lens.loc_s(:,2),mu_an,'-');
    %fLines = plot(lens.loc_s(:,2),(mu_an - mu_im)./mu_an,'o');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\mu(x_s)'), title('Fractional error in magnification factor');
end

% numerical plots
subplot(2,2,1); hold all;
plot(lens.loc_s(:,2),rho_im,'.');

subplot(2,2,2); hold all;
plot(lens.loc_s(:,2),tau_im,'.');

subplot(2,2,3); hold all;
plot(lens.loc_s(:,2),mu_im,'.');

% subplot(2,2,4); 
% plot(sourceloc(:,1),sourceloc(:,2),'ks');
% hold all;
% polar(phi_im,rho_im,'o');
% xlabel(lens.names.scale); ylabel(lens.names.scale); title('Source location');
%{
if(flag_plot==2||3)
    %subplot(nrows,ncols,1)
    %xlabel(lens.names.scale); ylabel(lens.names.scale);
    %subplot(nrows,ncols,2)
    %xlabel(lens.names.scale); ylabel('Scaled time delay');
    mtit(char(['Time delay surface: ', ...
        lens.names.title], lens.names.paramstr));
    filename = strcat('timesurf-',lens.names.filestr,'-new-all');
    filename = regexprep(filename,'+',''); % + forbidden in unix filenames
    saveas(gcf,filename,'png'); %saveas(gcf,filename,'fig');
end
%}
%% Functions to define y, dydx, d2ydx2

function y = y_fn(x)
X = abs(x);
% 2. mass distribution m(|x|)
switch(num_model)
    case{1}
        m = ones(size(x));
        
    case{2} % m/x = x/x_max^2; 1/x
        m(X < x_max) = x(X < x_max).^2./x_max.^2;
        m(X > x_max) = 1;
        
    case{3,4}
        g(X > 1) = 2./sqrt(X(X > 1).^2 - 1) .* atan(sqrt((X(X > 1) - 1)./(X(X > 1) + 1)));
        g(X < 1) = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
        g(X == 1) = 1;
        m = 4*kappa_s*(log(X/2) + g);
end
% rho_s = f(x) same as 0 = f(x) - rho_s
y = x - m./x - rho_s;
end % function

function dydx = dydx_fn(x)
X = abs(x);
switch(num_model)
    case{1}
        m = 1; dmdx = 0;
        
    case{2}
        m(X <= x_max) = x(X <= x_max).^2./x_max.^2;
        dmdx = 2.*x(X <= x_max)./x_max^2;
        m(X > x_max) = 1;
        dmdx(X > x_max) = 0;
        
    case{3,4}
        ind = find(X > 1);
        g(ind) = 2./sqrt(X(ind).^2-1) .* atan(sqrt((X(ind)-1)./(X(ind)+1)));
        dg(ind) = -2.*(sqrt((X(ind)-1)./(X(ind)+1)) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1))) + X(ind)-1) ./ ((X(ind)-1).^2.*X(ind).*(X(ind)+1));
        clear ind
        ind = find(X < 1);
        g(ind) = 2./sqrt(1-X(ind).^2) .* atanh(sqrt((1-X(ind))./(1+X(ind))));
        dg(ind) = 2.*(sqrt((1-X(ind))./(X(ind)+1)) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) + X(ind)-1) ./ (X(ind).*sqrt((1-X(ind))./(X(ind)+1)) .* (1-X(ind).^2).^(3/2));
        clear ind
        g(X == 1) = 1;
        dg(X == 1) = 0;
        
        m = 4*kappa_s*(log(X./2) + g);
        dmdx = 4*kappa_s*(1./(X./2) + dg);
end
dydx = 1 - (x.*dmdx - m)./x.^2;
end % function

function d2ydx2 = d2ydx2_fn(x)
X = abs(x);
switch(num_model)
    case{1}
        m = 1; dmdx = 0; d2mdx2 = 0;
        
    case{2}
        m(X <= x_max) = x(X <= x_max).^2./x_max.^2;
        dmdx = 2.*x(X <= x_max)./x_max^2;
        d2mdx2 = 2/x_max^2;
        m(X > x_max) = 1;
        dmdx(X > x_max) = 0;
        d2mdx2(X > x_max) = 0;
        
    case{3,4}
        ind = find(X > 1);
        g(ind) = 2./sqrt(X(ind).^2-1) .* atan(sqrt((X(ind)-1)./(X(ind)+1)));
        dg(ind) = -2*(sqrt((X(ind)-1)./(X(ind)+1)) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1)))+X(ind)-1) ./ ((X(ind)-1).^2.*X(ind).*(X(ind)+1));
        d2g(ind) = (-4*X(ind).^3 + 4*X(ind).^2 + 2*sqrt((X(ind)-1)./(X(ind)+1)) .* (2*X(ind).^2+1) .* X(ind).^2 .* atan(sqrt((X(ind)-1)./(X(ind)+1))) +X(ind)-1) ./ (X(ind).^2 .* sqrt((X(ind)-1)./(X(ind)+1)) .* (X(ind).^2-1).^(5/2));
        clear ind
        ind = find(X < 1);
        g(ind) = 2./sqrt(1-X(ind).^2) .* atanh(sqrt((1-X(ind))./(1+X(ind))));
        dg(ind) = 2.*(sqrt((1-X(ind))./(X(ind)+1)) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) + X(ind)-1) ./ (X(ind).*sqrt((1-X(ind))./(X(ind) + 1)) .* (1-X(ind).^2).^(3/2));
        d2g(ind) = (4.*X(ind).^3 - 4.*X(ind).^2 + 2.*sqrt((1-X(ind))./(X(ind)+1)) .* (2*X(ind).^2+1) .* X(ind).^2 .* atanh(sqrt((1-X(ind))./(X(ind)+1))) - X(ind)+1) ./ (X(ind).^2 .* sqrt((1-X(ind))./(X(ind)+1)) .* (1-X(ind).^2).^(5/2));
        clear ind
        g(X == 1) = 1;
        dg(X == 1) = 0;
        d2g(X == 1) = 0;
        
        m = 4*kappa_s*(log(X./2) + g);
        dmdx = 4*kappa_s*(2./X + dg);
        d2mdx2 = 4*kappa_s*(-2./X.^2 + d2g);
end
d2ydx2 = (-2.*m + 2.*x.*dmdx - d2mdx2)./x.^3;
end % function

function mu_x = mu_fn(x)
X = abs(x); m = ones(size(x)); kappa_x = zeros(size(x));

switch(num_model)
    case{1}
        m = 1; kappa_x(ceil(length(x)/2)) = 1;
        
    case{2} % m/x = x/x_max^2; 1/x
        m(X < x_max) = x(X < x_max).^2./x_max.^2;
        kappa_x(X <= x_max) = kappa_0 * pix^2; 
        
    case{3,4}
        g=zeros(size(x));
        g(X > 1) = 2./sqrt(X(X > 1).^2 - 1) .* atan(sqrt((X(X > 1) - 1)./(X(X > 1) + 1)));
        g(X < 1) = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
        g(X == 1) = 1;
        m = 4*kappa_s*(log(X/2) + g);
        
        f=zeros(size(x));
        f(X > 1) = 1 - 2./sqrt(X(X > 1).^2 - 1).*atan(sqrt((X(X > 1) - 1)./(1 + X(X > 1))));
        f(X < 1) = 1 - 2./sqrt(1 - X(X < 1).^2).*atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
        f(X ==1) = 0;
        kappa_x = 2*kappa_s.*f./(X.^2 - 1);
end
gamma_x = m./x.^2 - kappa_x; % shear: eqn. 8.15 SE&F; axi-symmetric only
detA = (1 - kappa_x).^2 - gamma_x.^2; % eqn. 5.25 SE&F; general
mu_x = abs(detA.^(-1)); % |inverse| = magnification
end % function

end % function