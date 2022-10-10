function plot_summary(lens,dist_mode)
%% Plot brightness and time delay curves during transit
global pc

Dd = dist_mode{2}; Ds = dist_mode{3};
len = size(lens.loc_s,1);

tau_im = lens.tau_im;
mu_im  = lens.mu_im;
x_im   = lens.x_im;
phi_im = lens.phi_im;
x_s    = lens.loc_s(:,2);
[sourceloc(:,1) sourceloc(:,2)] = pol2cart(lens.loc_s(:,1),lens.loc_s(:,2));

% convert from cylindrical polars to "lensing" polars
ix1 = find(phi_im ~= repmat(lens.loc_s(:,1),[1 size(phi_im,2)]));
ix2 = find(lens.loc_s(:,1) >= pi/2 | lens.loc_s(:,1) < -pi/2);
x_im(ix1) = -x_im(ix1); x_s(ix2) = -x_s(ix2); clear ix* % set r = -r

% colour
map = [0 0 0; 0.5 0 0.5; 0 0 0.9; 0 1 1; 0 1 0; 1 1 0; 1 0 0]; % hsv2
map = colormap_helper(map,len);
colour = map; colour = brighten(colour,-.25);
colour_cell = num2cell(colour,2);

% plotting
figure('visible','on');

% analytical plots, if soln extant
if(isfield(lens,'rho_an'))
    rho_an = cellextract(lens.rho_an);
    subplot(2,2,1); hold all;
    fLines = plot(x_s,rho_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('r(x_s)'),
    title('Fractional error in radial image locations');
end; hold off

if(isfield(lens,'tau_an'))
    tau_an = cellextract(lens.tau_an);
    subplot(2,2,2); hold all;
    fLines = plot(x_s,tau_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Fractional error in time delay');
end; hold off

if(isfield(lens,'mu_an'))
    mu_an = cellextract(lens.mu_an);
    subplot(2,2,3); hold all;
    fLines = plot(x_s,mu_an,'-');
    fGroup = hggroup; set(fLines,'Parent',fGroup);
    set(get(get(fGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(fGroup,'DisplayName','analytical soln.');
    xlabel(lens.names.scale), ylabel('\mu(x_s)'), title('Fractional error in magnification factor');
end; hold off

% numeric plots
subplot(2,2,1), hold all;
fHandles = plot(x_s,x_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('y(x_s)'), title('Radial image locations');
legend('Location','best'), axis tight; hold off;

subplot(2,2,2), hold all;
fHandles = plot(x_s,tau_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('\tau(x_s)'), title('Absolute time delay');
legend('Location','best'), axis tight; hold off;

subplot(2,2,3), hold all;
fHandles = plot(x_s,mu_im,'.');
fGroup = hggroup; set(fHandles,'Parent',fGroup);
set(get(get(fGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(fGroup,'DisplayName','numerical soln.');
xlabel(lens.names.scale), ylabel('\mu(x_s)'), title('Magnification factor');
legend('Location','best'), axis tight; hold off;

subplot(2,2,4); hold all;
scatter(sourceloc(:,1),sourceloc(:,2),36,colour,'Marker','s','DisplayName','lens');
parfor i=1:len
    hCLines(i,:) = polar(phi_im(i,:),abs(x_im(i,:)),'.');
end; clear i % for
hCGroup = hggroup; set(hCLines,{'Color'}, colour_cell,'Parent',hCGroup);
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
set(hCGroup,'DisplayName','images');
xlabel(lens.names.scale); ylabel(lens.names.scale); title('Location in lens plane');
axis tight; legend('Location','best');
hold off; clear h*

% Title for all plots incl. params and dist
title_str(1) = {['Source transiting ',lens.names.title]};
title_str(2) = {lens.names.paramstr};
title_str(3) = {['D_d = ',num2str(Dd/pc,'%8.2g'),' pc, D_s = ',num2str(Ds/pc,'%8.2g'),' pc']};

%mtit(char(title_str{1:2}));

%{
    if(flag_plot==2||3)
        filename = strcat('summary-',lens.names.filestr);
        filename = regexprep(filename,'+',''); % + forbidden in unix filenames
        saveas(gcf,filename,'png'); %saveas(gcf,filename,'fig');
    end %if
%}
end % function
