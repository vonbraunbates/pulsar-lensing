%{
% multiple_lenses.m fn(dist_cell,M_arr,num_model,zeta_arr,y_arr)

% INPUT:
% distance to lens, source plane(s) dist_cell {z1,...,zk,zs} or Di if z ~ 0
% masses of each lens M_arr = [M1 ... Mk] in M_sol or M_earth
% num_model (same for ALL lenses) 0, 1, 2
% lens locations x(l,1:2) is lth lens relative to that lens plane
% source locations y(k,1:2) is kth source rel. tp optical axis
%
% OUTPUT:
% image locations x1(y) = radial co-ords of images (r < 0 --> phi = -phi_s)
% magnification mu(x1;y): total, calculated recursively
% time delays tau(x1;y): total, scaled
% plot of pulsar signal T.*mu NOT scaled (in s)
%
%}
function result = multiple_lenses(dist_arr,M_arr,num_model,x_l,y_arr)
global c G pc t_H
% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
t_H = 3.09E17/0.6; % ! h=0.6 in Schneider

% 0. Setup
N   = length(M_arr); % nr. of lens planes
len = size(y_arr,1); % one colour per source location
map = [0 0 0; 0.5 0 0.5; 0 0 0.9; 0 1 1; 0 1 0; 1 1 0; 1 0 0]; % hsv2
map = colormap_helper(map,len);
colour = map; colour = brighten(colour,-.25);

% 1. Scaling
[z D] = tdelay(dist_arr);
% Only need Ds = D(1,k+1), Dd = D(1,k), Dds = D(k,k+1)
Ds = D(1,2:N+1)'; Dd = D(1,1:N)'; Dds = diag(D); Dds = Dds(2:end);
M_mks = M_arr'.*M_earth; %3e-2*M_sol;          % total lens mass
r_E = sqrt(4*G) .* sqrt(M_mks./c^2)' .* sqrt(Dd.*Dds./Ds); % Einstein radii

% 3. Lens parameters dependsing on model
switch(num_model)
    case{1}
        params = [M_arr(:), r_E(:)];
        
    case{2}
        if(~scale_mode)
            r_max = 10^(-4.5)*pc * ones([1 N_lens]); % lens radius
            x_max = r_max./r_E;
        else
            x_max = .1*ones([1 N_lens]);                          % r_max/zeta_0
        end % if
        params = [M_arr(:), r_E(:), x_max(:)];
        
    case{3} % NFW
        r_0 = 10^(-4.5)*pc*ones([1 N_lens]); % scale radius of lens (pc)
        if(~scale_mode)
            sigma_cr = c^2./(4*pi*G) .* Ds./(Dd.*Dds); % kg/m^2
            r_max = 10^(-3)*pc; c_nfw = r_max./r_0; x_max = c_nfw; 
            delta_c = M_arr./((4*pi*r_0.^3).*(log(1 + c_nfw) - c_nfw./(1+c_nfw)))
            kappa_s = delta_c.*r_s./sigma_cr
        end
        
        params = [M_arr(:), r_0(:), x_max(:), kappa_s(:)];
        clear rho_crit r_s
end % switch

% 2. Recursive mu, tau for each lens plane
% For each original lens position, we require a row vector of images in the
% final lens plane with tau(x1;y) mu(x1;y) rho(x1;y)
[temp(:,1) temp(:,2)] = cart2pol(y_arr(:,1), y_arr(:,2));
result.loc_s = temp; clear temp;
for l = 1:size(y_arr,1)
     source = y_arr(l,:); x_im = source; j = N+1; x_scale = 1; 
% x_im = y_arr; j = N+1; x_scale = 1; l=1; y_ind = [1 2]';
    % Trace each initial source through the lens planes
    for k = N:-1:1
        % k-th plane is lens plane, j-th plane is source plane, where
        % j = min{k+1, ... , N: M(j)>0}, i.e. nearest plane with mass
        dist_mode = {dist_arr{1},dist_arr{k+1},dist_arr{j+1}};
        if(M_arr(k) > 0)
            % co-ord transformation [x^(j) := y^(k)]
            if(j < N+1) % if j-th plane is not S
                zeta_j  = params(j,2); % x^(j) = zeta^(j)/zeta_0
                eta_k   = params(k,2)*Dd(j)/Dd(k); % eta_0 = zeta_0*Ds/Dd
                x_scale = zeta_j/eta_k % y^(k) = eta^(k)/eta_0
            end % if
            % x_s = x_im[x^(j) := y^(k) co-ords] + [shift rel to new axis]
            source = x_im ./ x_scale + repmat(-x_l(k,:),[size(x_im,1) 1]);
            lens = newlensinght(source,dist_mode,params(k,:),num_model,0,1,1,1);
            % get 1 cell entry per lens field
            tau_im   = lens.tau_im';
            mu_im{k} = lens.mu_im';
            rho_im   = lens.x_im';
            phi_im   = lens.phi_im';
%             y_ind    = repmat(y_ind',[size(rho_im,1) 1]);
%             y_ind    = y_ind(:);
            % Translate y=abs(x_s) relative to optical axis (x_s in R, y in R+)
            rho_s    = hypot(source(:,1) + x_l(k,1), source(:,2) + x_l(k,2));
            % Total time delay relative to optical axis
            tau_opt  = tau_im + repmat(0.5*rho_s.^2,[1 size(tau_im,2)]);
            % scale time delays between each lens plane pair
            % T(x,y) eq. 5.45 Schneider (zeta_0 = params(k,2))
            t_opt{k} = params(k,2)^2/c * Ds(k)/(Dd(k)*Dds(k)) * (1+z(k)).*tau_opt;
            % Translate x_im relative to opt. axis
            clear rho_s tau_opt x_im % resets size of x_im to 0x0
            [x_im(:,1), x_im(:,2)] = pol2cart(phi_im(:),rho_im(:));
            x_im = x_im + repmat(x_l(k,:),[size(x_im,1) 1]);
            % This will be a source plane for the next lens plane
            j = k;
        end
    end; clear k % for
    result.names.title    = lens.names.title;
    result.names.paramstr = lens.names.paramstr;
    result.names.scale    = lens.names.scale;
    % clear lens
    % Images in final lens plane, Cartesian co-ords
    result.x_im{l} = x_im;
    % Total time delay is sum over all lens planes -- for some reason this
    % result is already scaled!
    % Also have to propagate delys between planes (array sizes change)
    result.tau_im{l} = sum(cellextract(t_opt,[],0))';
    % This result is probably shifty too, if the tds are wrong?
    result.mu_im{l}  = sum(cellextract(mu_im,[],0))';
    % Plot corresponding summary in scaled units
    % This does something odd to the image locations!!
    plot_summary(result,dist_mode);
    % Title for all plots incl. params and dist
    title_str(1) = {['Source transiting ',result.names.title]};
    title_str(2) = {result.names.paramstr};
    title_str(3) = {['Lens distances: ',sprintf('%8.2g',dist_arr{2:N+1}),'kpc. Source distance:',sprintf('%8.2g',dist_arr{N+2})]};
end; clear l % for

% 2. Effect of time delay on pulsar signal - this is broken!
%%{
% create new signal to account for lensing
t_delay = cellextract(result.tau_im);  clear t;
mu_tot  = cellextract(result.mu_im); clear mu;
siz     = size(t_delay,2);                  % largest nr of images
period  = 1e-6;                             % microsecond pulsar
t_old   = linspace(0, period*len, len)';    % pulse at each sourceloc
t_new   = repmat(t_old,[1 siz]) + t_delay;  % replicate to nr. of images
s_old   = ones([len 1]);                    % normalised amplitude
s_new   = repmat(s_old,[1 siz]) .* mu_tot;  % signal amplitude
% plot
figure('visible','on'); hold on;
hSLines = scatter(t_old,s_old,36,colour,'Marker','*','DisplayName','source');
parfor i=1:len
    hILines(i,:) = plot(t_new(i,:),s_new(i,:),'o:','Color',colour(i,:));
end; clear i % for
hCGroup = hggroup; set(hILines,'Parent',hCGroup);
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
set(hCGroup,'DisplayName','images per source location');
xlabel('Time (s)'); ylabel('Flux relative to source'); legend('Location','Best');
axis tight; mtit(char(title_str)); % for A \n B C use char(A,[B,C]);
%}
end % function