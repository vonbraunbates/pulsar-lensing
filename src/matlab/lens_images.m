function lens = lens_images(sourceloc,num_model,params,varargin)

global c G pc t_H eps_ % mks universal constants; single precision machine eps

if(nargin==3) 
    flag_plot    = 0;  % plot time delay surf, cross-sections
    flag_summary = 0;  % plot scaled tau(x), mu(x), rho(x) as source moves
    scale_mode   = 1;  % if 1, use scaled lens params; if 0 use mks units
elseif(nargin==4)
    flag_plot    = varargin{1}; 
    if(flag_plot); x1 = varargin{3}; x2 = varargin{4}; end
    flag_summary = varargin{2}; 
end % if
allcolors = get(0,'defaultAxesColorOrder');
allcolors(1,3) = 0.75; allcolors(3,1) = 0.75; % change 'red', 'blue'

%% subplot setup
len = size(sourceloc,1); % nr of plots
if(flag_plot)
    if(flag_plot==1||len==2)
        nrows=1; ncols=3;  % 3 subplots
    elseif(flag_plot==2)
        nrows = ceil(len/4); ncols = 4;
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
    
    clear scrsz ax.min ax.max ax.size ax.box
end

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

%% Lens and model parameters


%% Set up lens plane and convolution
num = pow2(8); % num pixels in [0,max_radius] 
max_radius = 25; % box edge 
pix = max_radius / num  % pixel size in scaled units
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

sigma_fft   = ht(kappa);
if(num_model==1) % kappa(r) = dirac(r)
    sigma_fft = 1/(2*pi)*ones(num,1); % manual HT
end % if
krn_fft     = ht(kernel);
psi_fft     = 2*pi*sigma_fft .* krn_fft; % H[f1**f2] = 2pi.F1.F2
psi         = real ( iht( psi_fft ) ) ;     % psi = 1/pi f**g
%%{ 
% rotate to get 2d potential time delay
if(flag_plot);
    [polyrev.X,polyrev.Y,polyrev.Z] = polyrevolve(rho,psi(1:num),0.4);
    polyrev.F = TriScatteredInterp(polyrev.X,polyrev.Y,polyrev.Z);
    td_pot    = polyrev.F(x1,x2);
end;
%}
for i = 1:len
    source = sourceloc(i,:); % dimensionless source loc.
    [phi_s rho_s] = cart2pol(source(1),source(2));
    sourcestr = sprintf('%0.5e %0.5e',source);
    
    % geometric time delay
    if(flag_plot); td_geom = 0.5.*((x1-source(1)).^2 + (x2-source(2)).^2); end;
    %[xi,yi] = pol2cart(phi_s*ones(size(x)),x); % grid (phi,rho) --> (x,y)
    %lens.tau{i}(:,3) = 0.5.*((xi-source(1)).^2 + (yi-source(2)).^2);
    lens.tau{i}(:,3) = 0.5.*(x - rho_s).^2;
    
    % total time delay
    lens.tau{i}(:,2) = [flipud(psi(:));psi(:)]; 
    if(flag_plot); td_tot = td_geom - td_pot; end;
    lens.tau{i}(:,1) = lens.tau{i}(:,3) - lens.tau{i}(:,2);
    % magnification factor
    lens.mu{i} = mu_fn(x);
    
    % image locations:
    % 1. mass distribution m(|x|)
    % This goes in subfunctions at end:
    
    % 2. Solve lens equation NUMERICALLY
    x0 = rootsearch(@y_fn,@dydx_fn,@d2ydx2_fn,min(x),max(x))';
    
    % 3. polar co-ords
    lens.loc_im{i}(1,:) = x0;%sort(x0); 
    lens.loc_im{i}(2,:) = phi_s*ones(size(x0));
    
    % 4. interpolated solutions at image locations
    lens.tau_im{i} = interp1(x,lens.tau{i}(:,1),x0,'spline')';
    lens.mu_im{i}  = interp1(x,lens.mu{i},x0,'spline')';
    
    %% Analytical soln.
    % radial co-ord rho; source pos. y
    y = rho_s;

    switch(num_model)
        case{1}
            phi_an = phi_s*ones(2,1); lens.rho_an{i} = zeros(2,1);
            x_im(1) = (y + sqrt( y^2 + 4 ) )/2 ;
            x_im(2) = (y - sqrt( y^2 + 4 ) )/2 ;
            %x_im = sort(x_im);
            lens.rho_an{i} = x_im;
            lens.tau_an{i} = 0.5*(x_im - y).^2 - log(abs(x_im)) + log(abs(max(x)));
            lens.mu_an{i}  = abs((1 - x_im.^(-4)).^(-1));
            clear *_im
            
        case{2}
            if(x_max < 1)
                if(y < (1 - x_max^2)/x_max)
                    disp('y < (1 - x_max^2)/x_max')
                    phi_an = phi_s*ones(3,1); lens.rho_an{i} = zeros(3,1);
                    x_im(1) = y * x_max^2/(x_max^2 - 1) ;
                    x_im(2) = (y + sqrt( y^2 + 4 ) )/2 ;
                    x_im(3) = (y - sqrt( y^2 + 4 ) )/2 ;
                elseif(y == (1 - x_max^2)/x_max)
                    phi_an = phi_s*ones(2,1); lens.rho_an{i} = zeros(2,1);
                    x_im(1) = x_max ;
                    x_im(2) = 1/x_max ;
                elseif(y > (1 - x_max^2)/x_max)
                    phi_an = phi_s;
                    x_im = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            elseif(x_max > 1)
                phi_an = phi_s;
                if(y <= (x_max^2 - 1)/x_max)
                    x_im = y * x_max^2/(x_max^2 - 1) ;
                else
                    x_im = (y + sqrt( y^2 + 4 ) )/2 ;
                end
            else % if x_max = 1 all points are focused onto the optical axis
                phi_an = phi_s; x_im = 0;
            end
            % SMD kappa = surface mass per pix = surface density * pix. area
            % So x_0^(-2) = kappa_0 --> kappa_0*pix^2
            mu_im = ((1 - 1./x_im.^2).*(1 + 1./x_im.^2)).^(-1); % true x_im > x0 
            [lens.rho_an{i},ix] = sort(x_im); 
            lens.mu_an{i} = abs(mu_im(ix));
            ind = find(abs(lens.rho_an{i}) < x_max);
            lens.mu_an{i}(ind) = (1 - kappa_0).^(-2); % for x_im < x0
            clear *_im
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
        for l=1:size(lens.loc_im{i},2)
            hSLines(l) = polar(lens.loc_im{i}(2,l),lens.loc_im{i}(1,l));
            set(hSLines(l),'Marker','x','color',allcolors(l,:));
        end; clear l
        hSGroup = hggroup; set(hSLines,'Parent',hSGroup);
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(hSGroup,'DisplayName','Numerical image location');
        
        % plot analytical (o) image locations
        if(isfield(lens,'rho_an'))
            for l=1:size(lens.rho_an{i},2)
                hCLines(l) = polar(phi_an(l),lens.rho_an{i}(l));
                set(hCLines(l),'Marker','o','color',allcolors(l,:));
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
        plot(x,lens.tau{i}(:,1),'k-', ...
            x,lens.tau{i}(:,2),'k-.', ...
            x,lens.tau{i}(:,3),'k:');
        for l=1:size(lens.tau_im{i},1)
            hSLines(l) = plot([-num/2*pix num/2*pix], ...
                [lens.tau_im{i}(l) lens.tau_im{i}(l)], '--','color',allcolors(l,:));
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
        plot(x,lens.mu{i},'k-','DisplayName','Magnification factor');
        
        % plot numerical (x) image magnification
        for l=1:size(lens.loc_im{i},2)
            hSLines(l) = plot(lens.loc_im{i}(1,l),lens.mu_im{i}(l),'x','color',allcolors(l,:));
             if(isfield(lens,'mu_an'))
                 fprintf('analytical value: %.5f ',lens.mu_an{i}(l));
             end
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
    end % if

end; clear i % for

%% Plot brightness and time delay curves during transit

% introduce negative radii
ix = find(lens.loc_s(:,1) < 0);             % for phi < 0
lens.loc_s(ix,2) = -lens.loc_s(ix,2);       % set r = -r
lens.loc_s(ix,1) = lens.loc_s(ix,1) - pi;   % set phi = phi - pi
clear ix

% plot does not work with cell types 
tau_im = cellextract(lens.tau_im);
mu_im  = cellextract(lens.mu_im);
rho_im = cellextract(lens.loc_im,'@(x)x(1,:)');
phi_im = cellextract(lens.loc_im,'@(x)x(2,:)');
y = abs(lens.loc_s(:,2));

% colour
map = [0 0 0; 0.5 0 0.5; 0 0 0.9; 0 1 1; 0 1 0; 1 1 0; 1 0 0]; % hsv2
map = colormap_helper(map,len);
colour = map; colour = brighten(colour,-.25);
colour_cell = num2cell(colour,2);
%%{
figure('visible','on');%figure(4);

% analytical plots, if soln extant
if(isfield(lens,'rho_an'))
    rho_an = cellextract(lens.rho_an);
    subplot(2,2,1); hold all;
    fLines = plot(lens.loc_s(:,2),rho_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('r(x_s)'), title('Fractional error in radial image locations');
end; hold off

if(isfield(lens,'tau_an'))
    tau_an = cellextract(lens.tau_an);
    subplot(2,2,2); hold all;
    fLines = plot(lens.loc_s(:,2),tau_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Fractional error in time delay');
end; hold off

if(isfield(lens,'mu_an'))
    mu_an = cellextract(lens.mu_an);
    subplot(2,2,3); hold all;
    fLines = plot(lens.loc_s(:,2),mu_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\mu(x_s)'), title('Fractional error in magnification factor');
end; hold off

% numeric plots

subplot(2,2,1), hold all;
fHandles = plot(lens.loc_s(:,2),rho_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('y(x_s)'), title('Radial image locations');
legend('Location','best'), axis tight; hold off;

subplot(2,2,2), hold all;
fHandles = plot(lens.loc_s(:,2),tau_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Absolute time delay');
legend('Location','best'), axis tight; hold off;

subplot(2,2,3), hold all;
fHandles = plot(lens.loc_s(:,2),mu_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('\mu(x_s)'), title('Magnification factor');
legend('Location','best'), axis tight; hold off;

subplot(2,2,4); hold all;
hSLines = scatter(sourceloc(:,1),sourceloc(:,2),36,colour,'Marker','s','DisplayName','source');
for i=1:len
    hCLines(i,:) = polar(phi_im(i,:),rho_im(i,:),'.');
end; clear i % for
hCGroup = hggroup; set(hCLines,{'Color'}, colour_cell,'Parent',hCGroup);
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
set(hCGroup,'DisplayName','images');
xlabel(lens.names.scale); ylabel(lens.names.scale); title('Location in lens plane');
axis tight; legend('Location','best');
hold off; clear h*
%}
% Title for all plots
title_str(1) = {['Source transiting ',lens.names.title]};
title_str(2) = {[lens.names.paramstr,'N = ',num2str(num),', R = ',num2str(max_radius)]};
mtit(char(title_str));

%{
if(flag_plot==2||3)
    filename = strcat('summary-',lens.names.filestr,'-new-all');
    filename = regexprep(filename,'+',''); % + forbidden in unix filenames
    saveas(gcf,filename,'png'); %saveas(gcf,filename,'fig');
end %if
%}
%}

%% Impose time delay on signal
% Total time delay relative to optical axis
tau = tau_im + repmat(0.5*y.^2,[1 size(tau_im,2)]);

% Scaling
t = tdelay(dist_mode,tau,lens.loc_s(:,2),zeta_0);


%% Surface mass density
    function kappa = kappa_fn(rho,num_model)
        kappa = zeros(size(rho)); % convergence (SMD) without zp (at this stage!)
        
        switch(num_model)
            case{1}
%                 % \int_R dx dirac(x) = 1 iff dirac(0) = 1/(\int_R dx)
%                 % normally R is the central pixel, area 1/(pix)^2 BUT
%                 % we want independence from grid fineness so set area(R)=1
%                 dirac_x = 1; 
%                 % dirac(at) = |a|^(-n) * dirac(t) in R^n
%                 dirac_zeta = (zeta_0^(-2)) * dirac_x;
%                 % change coords: dirac(x,y) = dirac(rho)/(2*pi*rho)
%                 dirac_zeta(num) = 0;
%                 dirac_zeta = 1./rho' .* dirac_zeta./(2*pi);
%                 % kappa = sigma(zeta)/sigma_cr = M*dirac(zeta)/sigma_cr
%                 kappa = dirac_zeta/sigma_cr;
                  dirac_r(1)=pi; dirac_r(num)=0; % Cartesian co-ords
                  kappa = M*dirac_r(:)./(2*pi*rho);   % Cylindrical co-ords
                
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
        
       % kappa(2*num) = 0; % zero-padding
    end % function

%% Magnification factor
    function mu_x = mu_fn(x)
        X = abs(x); m = ones(size(x)); kappa_x = zeros(size(x));
        
        switch(num_model)
            case{1}
                m = ones(size(x)); 
                % zero except at x=0
                kappa_x(ceil(length(x)/2)) = 1;
                
            case{2} % m/x = x/x_max^2; 1/x
                m(X < x_max) = x(X < x_max).^2./x_max.^2;
                kappa_x(X <= x_max) = x_max^(-2);%kappa_0 * pix^2;
                
            case{3,4}
                g=zeros(size(x));
                g(X > 1) = 2./sqrt(X(X > 1).^2 - 1) .* atan(sqrt((X(X > 1) - 1)./(X(X > 1) + 1)));
                g(X < 1) = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                g(X ==1) = 1;
                m = 4*kappa_s*(log(X/2) + g);
                
                f=zeros(size(x));
                f(X > 1) = 1 - 2./sqrt(X(X > 1).^2 - 1).*atan(sqrt((X(X > 1) - 1)./(1 + X(X > 1))));
                f(X < 1) = 1 - 2./sqrt(1 - X(X < 1).^2).*atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                f(X ==1) = 0;
                kappa_x = 2*kappa_s.*f./(X.^2 - 1);
        end
        gamma_x = m./(x.^2) - kappa_x; % shear: eqn. 8.15 SE&F; axi-symmetric only
        detA1 = (1 - m./x.^2).*(1 + m./x.^2 - 2*kappa_x); % eqn. 8.16 SE&F; 
        %detA2 = (1 - kappa_x).^2 - gamma_x.^2; % eqn. 5.25 SE&F; general
        %if(abs(detA1 - detA2) < eps_)
            mu_x = abs(detA1.^(-1)); % |inverse| = magnification
        %else
        %    disp('Error: detA formulae do not agree.')
        %end; % if
    end % function

%% Functions to define y, dydx, d2ydx2

    function y = y_fn(x)
        X = abs(x);
        % 2. mass distribution m(|x|)
        switch(num_model)
            case{1}
                m = ones(size(x));
                
            case{2} % m/x = x/x_max^2; 1/x
                m(X <= x_max) = x(X <= x_max).^2./x_max.^2;
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
                m = ones(size(x)); dmdx = zeros(size(x));
                
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
                m = ones(size(x)); dmdx = zeros(size(x)); d2mdx2 = zeros(size(x));
                
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

end % function