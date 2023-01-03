function out_struct = lensingdb(dist_cell,M_arr,num_model,scale_mode,varargin)

%{
lensingdb.m
Catalogue of lensing signals for a given set of parameter space
INPUT:
num_model = lensing model
M_arr = array of mass lenses (in M_earth units)
dist_arr = {'mode',d1,...dn} where 'mode' is either 'z' or 'dist'
params = cell of parameters for each model
OUTPUT:
x   = radial position array of all images, lensing polar co-ords
phi = angular position array of all images, lensing polar co-ords
mu  = magnification factors for each image
tau = time delay relative to optical axis for each image
%}

% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
t_H = 3.09E17/0.72; % ! h=0.6 in Schneider

% 1. set variables
T_obs   = 3e7; % period between first, last observation (s); 1yr = 3e7s
N_pulse = 100; % number of signals in traversal of source plane
N_lens  = numel(M_arr); % number of lenses 
% impact parameter from normal dist in [0,1]
if(nargin < 5)
    b_arr = rand(1,N_lens); % assign b at random for varying Dd
elseif(nargin == 5 && ndims(varargin{:})==2); % v{:} is list of b
    b_arr = sort(varargin{:},'ascend');
    b_max = b_arr(end); b_arr = b_arr(1:end-1); 
elseif(nargin == 5 && ndims(varargin{:})==3); % v{:} is array of y
    y = varargin{:}; 
    b_arr = min(hypot(y(:,1,:),y(:,2,:)));
end % if

% 2. Scaling
Dd = reshape(dist_cell{2}*pc,[N_lens 1]); Ds = dist_cell{3}*pc; Dds = Ds - Dd;
z = [Dd; Ds]./(c*t_H); 
M_arr = reshape(M_arr.*M_earth,[N_lens 1]); % lens masses (kg)
r_E = sqrt(4*G) .* sqrt(M_arr./c^2) .* sqrt(Dd.*Dds./Ds); % Einstein radii

%% 3. Lens parameters depending on model
switch(num_model)
    case{1}
        params = [M_arr(:), r_E(:)];
        
    case{2}
        if(~scale_mode)
            r_max = 10^(-4.5)*pc * ones([1 N_lens]); % lens radius
            x_max = r_max./r_E;
        else
            x_max = .1*ones([1 N_lens]);             % r_max/zeta_0
        end % if
        params = [M_arr(:), r_E(:), x_max(:)];
        
    case{3} % NFW
        r_0 = 10^(-4.5)*pc*ones([N_lens 1]); % scale radius of lens (pc)
        if(~scale_mode)
            sigma_cr = c^2./(4*pi*G) .* Ds./(Dd.*Dds); % kg/m^2
            r_max = 10^(-3)*pc; c_nfw = r_max./r_0; x_max = c_nfw; 
            delta_c = M_arr./((4*pi*r_0.^3).*(log(1 + c_nfw) - c_nfw./(1+c_nfw)));
            kappa_s = delta_c.*r_0./sigma_cr;
        else
            x_max = 10*ones([1 N_lens]);  % extent of lens (>1 b/c units of r_s)
            kappa_s = .5*ones([1 N_lens]); % ~ 1 (Bartelmann)
        end
        
        params = [M_arr(:), r_0(:), x_max(:), kappa_s(:)];
        clear rho_crit r_s
end % switch

%% 4. Generate lens paths
v = 200*ones([N_lens 1]); %(v_mu - v_sigma) + 2*v_sigma.*rand(N_lens,1); % draw random velocities from dist.
% change this so T starts randomly!
T = linspace(-0.5*T_obs, 0.5*T_obs, N_pulse)';
if(~exist('y','var'))
% get vt line of transit (centred on origin)
y_len = T*(1e3*v'); 
eta_0 = repmat(params(:,2)'.*Ds./Dd(:)', [N_pulse 1]); % eta_0 = zeta_0*Ds/Dd
y_len = y_len./eta_0; % scaling
y_len = reshape(y_len,[N_pulse 1 N_lens]);
% impact parameter (closest approach to lens) is x_lens
b = reshape(repmat(b_arr,[2*N_pulse 1]), [N_pulse 2 N_lens]);
% pick random inclination for line [-pi,pi]
y_ang1 = repmat(2*pi*rand([1 1 N_lens]) - pi, [N_pulse 1 1]); % zeros([N_pulse 1 N_lens]);
% pick random rotation of line about x-axis
y_ang2 = repmat(2*pi*rand([1 1 N_lens]), [N_pulse 1 1]);
% polar to Cartesian co-ords;
y(:,1,:) = y_len.*cos(y_ang1) + b(:,1,:).*cos(y_ang2);
y(:,2,:) = y_len.*sin(y_ang1) + b(:,2,:).*sin(y_ang2);
%figure; plot(squeeze(y(:,1,:)),squeeze(y(:,2,:)),'.');
clear b v y_ang y_len
end % if

% Exclude sources outside box
mask(:,1,:) = (abs(y(:,1,:)) < b_max);
mask(:,2,:) = (abs(y(:,2,:)) < b_max);
mask(:,3,:) = ((mask(:,1,:) ==1) & (mask(:,2,:)==1));
    
%% 5. For each lens, find [mu tau]
for k = 1:N_lens
    dist_mode = {dist_cell{1},dist_cell{2}(k),dist_cell{3}};
    source_cell{k} = squeeze(y(mask(:,3,k),1:2,k)); 
    % shift co-ords so that lens is always at center
    lens   = newlensinght(-source_cell{k},dist_mode,params(k,:),num_model,0,0,0,1);
    tau_im = lens.tau_im;
    mu{k}  = lens.mu_im;
    % Translate y=abs(x_s) relative to optical axis (x_s in R, y in R+)
    rho_s  = hypot(source_cell{k}(:,1), source_cell{k}(:,2));
    % Total time delay relative to optical axis
    tau_opt = tau_im + repmat(0.5*rho_s.^2,[1 size(tau_im,2)]);
    % scale time delays between each lens plane pair
    % T(x,y) eq. 5.45 Schneider (zeta_0 = params(k,2))
    t{k} = params(k,2)^2/c * Ds/(Dd(k)*Dds(k)) * (1+z(k)).*tau_opt;
end; % for

%% 6. combine signals from different lenses and plot
t_struct = plot_tdelay(t,mu,T);
t_sorted = t_struct.t; s_sorted = t_struct.mu; T_res = t_struct.T_res; 
handles.fig = t_struct.fig; handles.ax = t_struct.ax;

% plot lens locations
set(handles.fig,'CurrentAxes',handles.ax(3));
plot(0,0,'ks','DisplayName','source'); hold on;
for j=1:N_lens
    scatter(source_cell{j}(:,1).*(eta_0(:,j)/pc), ...
        source_cell{j}(:,2).*(eta_0(:,j)/pc), 36, ...
        t_struct.colour{j},'Marker','.');
end; clear j % for
xlabel('pc','fontsize',10,'fontname','Century Schoolbook L'); 
ylabel('pc','fontsize',10,'fontname','Century Schoolbook L');  
title('Lens location (source plane)','fontsize',10,'fontname','Century Schoolbook L','fontweight','bold'); 

%% colourbar
if(isscalar(unique(dist_cell{2}))) % if lenses are all at same dist, cb is impact parameter
    b_unique = unique(b_arr); % sort
    if(N_lens == 1)
        legend('source',['lens b = ',num2str(b_arr,4)],'Location','Best'); % legend
    elseif(N_lens <= 10)
        cbticks = cellfun(@(x)(num2str(x,'%-5.3f')),num2cell(b_unique),'UniformOutput',0);
        num_cbt = numel(cbticks);
        colormap(t_struct.cbar_map); % must go before setting caxis
        caxis([1 num_cbt]); % must go before cbfit
        handles.cb = colorbar('YTickLabels',cbticks,'YLim',[1 num_cbt], ...
            'Location','WestOutside','fontsize',10,'fontname','Century Schoolbook L');
        set(handles.cb,'YTickMode','manual'); % must go after colourbar()
        set(handles.cb,'YTick',[1 + .5*(num_cbt-1)/num_cbt : (num_cbt-1)/num_cbt : num_cbt]);
        cblabel(handles.cb,'impact parameter','fontsize',10,'fontname','Century Schoolbook L');
    else
        ind = mod([1:length(b_unique)],5); % mod(x,y) get every yth tick in x
        cbticks = cellfun(@(x)(num2str(x,'%-5.3f')),num2cell(b_unique(ind==1)),'UniformOutput',0);
        num_cbt = numel(cbticks);
        colormap(t_struct.cbar_map); % must go before setting caxis
        caxis([1 num_cbt]); % must go before cbfit
        handles.cb = colorbar('YTickLabels',cbticks,'YLim',[1 num_cbt], ...
            'Location','WestOutside','fontsize',10,'fontname','Century Schoolbook L');
        set(handles.cb,'YTickMode','manual'); % must go after colourbar()
        set(handles.cb,'YTick',[1 + .5*(num_cbt-1)/num_cbt : (num_cbt-1)/num_cbt : num_cbt]);
        cblabel(handles.cb,'impact parameter','fontsize',10,'fontname','Century Schoolbook L');
    end % if
else % if Dd varies, cb is distance as a fraction of source dist
    Dd_unique = Dd./Ds;
    if(N_lens == 1)
        legend('source',['lens b = ',num2str(b_arr,4),'r_0; Dd = ',num2str(Dd,4),' pc'],'Location','Best'); % legend
    elseif(numel(Dd_unique) <= 10)
        cbticks = cellfun(@(x)(num2str(x,'%-5.3f')),num2cell(Dd_unique),'UniformOutput',0);
        num_cbt = numel(cbticks);
        colormap(t_struct.cbar_map); % must go before setting caxis
        caxis([1 num_cbt]); % must go before cbfit
        handles.cb = colorbar('YTickLabels',cbticks,'YLim',[1 num_cbt], ...
            'Location','WestOutside','fontsize',10,'fontname','Century Schoolbook L');
        set(handles.cb,'YTickMode','manual');
        set(handles.cb,'YTick',[1 + .5*(num_cbt-1)/num_cbt : (num_cbt-1)/num_cbt : num_cbt]);
        cblabel(handles.cb,'lens dist. (units of source dist.)','fontsize',10,'fontname','Century Schoolbook L');
    else
        ind = mod([1:length(Dd_unique)],5); % mod(x,y) get every yth tick in x
        cbticks = cellfun(@(x)(num2str(x,'%-5.3f')),num2cell(Dd_unique(ind==1)),'UniformOutput',0);
        num_cbt = numel(cbticks);
        colormap(t_struct.cbar_map); % must go before setting caxis
        caxis([1 num_cbt]); % must go before cbfit
        handles.cb = colorbar('YTickLabels',cbticks,'YLim',[1 num_cbt], ...
            'Location','WestOutside','fontsize',10,'fontname','Century Schoolbook L');
        set(handles.cb,'YTickMode','manual'); % must go after colourbar()
        set(handles.cb,'YTick',[1 + .5*(num_cbt-1)/num_cbt : (num_cbt-1)/num_cbt : num_cbt]);
        cblabel(handles.cb,'lens dist. (units of source dist.)','fontsize',10,'fontname','Century Schoolbook L');
    end % if
end % if
axis tight;

% main title
if(isscalar(unique(Dd)))
    dist_str = ['D_d ',num2str(min(Dd)/(1e3*pc),'%-0.2g')];
else
    dist_str = [num2str(min(Dd)/(1e3*pc),'%-0.2g'),' < D_d < ',num2str(max(Dd)/(1e3*pc),'%-0.2g')];
end % if
title_str = {[lens.names.paramstr,' ',dist_str,' kpc D_s ', ...
    num2str(Ds/(1e3*pc),'%-0.2g'),' kpc T_{res} ',num2str(T_res),' s']};
title_str = regexprep(title_str,'e\+0(\w*)',' \\times 10^{$1}'); % e+0x --> 10^{x}
title_str = regexprep(title_str,'e-0(\w*)',' \\times 10^{-$1}'); % e-0x --> 10^{-x}
handles.mtit = mtit(handles.ax,title_str{:},'fontsize',12,'fontname','Century Schoolbook L','fontweight','bold'); 
set(handles.mtit.th,'Position',[0.5 1.035 9.16]);

%% output
wait = input('Have you finished adjusting figure? y/n: ','s');
if(wait=='y')
    out_fields = {'T','T_res','y','mu','t','t_sorted','s_sorted','handles'};
    for i=1:numel(out_fields)
        out_struct.(out_fields{i}) = eval(out_fields{i});
    end % for
else
    return
end % if
end % function