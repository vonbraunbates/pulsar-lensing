function t_struct = plot_tdelay(t,mu,t_source)

% get time delays and magnification
N_lens    = numel(t); % t is a cell: t{i} = dt(ith lens)
N_pulse   = size(t{1},1);
N_images  = max(cellfun(@(x)size(x,2), t));
tau_array = cellextract(t)';  % t_arr(:,i)  = t{ith lens}(:,:)
mu_array  = cellextract(mu)'; % mu_arr(:,i) = mu{ith lens}(:,:)

% setup signal
T_res    = 1e-5; t_struct.T_res = T_res; % timing residual of pulsar (s); typically 1us 
t_source = t_source + t_source(end); % [0,T] not [-T/2,T/2]
s_source = ones(size(t_source)); % original pulsar signal

%% Create signals
% add lensing effects to signal
t_images = repmat(t_source, [N_images N_lens]) + tau_array; 
t_images = reshape(t_images,[N_pulse N_lens*N_images]);
s_images = repmat(s_source, [N_images N_lens]) .* mu_array; 
s_images = reshape(s_images,[N_pulse N_lens*N_images]);
clear mu_array 

% bin signals by time
[temp1,ind] = sort(t_images(:)); % sort t
temp2 = s_images(:); temp2 = temp2(ind); % sort mu by t value
clear ind; 
k = 1; j = 1;
while k < length(temp1);
    % find signals close together
    dt = temp1 - temp1(k);
    t = temp1((dt >= 0) & (dt < 10*T_res));
    s = temp2((dt >= 0) & (dt < 10*T_res));
    k = find(dt > 10*T_res,1,'first');
    % bin only those signals
    [counts, bin] = histc(t,[min(t) : T_res : max(t)+T_res]);
    max_counts = max(counts); % largest nr of superposed signals
    l=1; % number of non-empty bins
    % Non-empty bins contain signals which will be superimposed
    for i = 1:length(counts)
        if(~isempty(t(bin==i)))
            t_mat(l,:) = vec2mat(t(bin==i),max_counts,NaN);
            s_mat(l,:) = vec2mat(s(bin==i),max_counts,0);
            l = l+1;
        end % if
    end; % for
    clear bin i max_counts nbins
    % get time, signal for composite
    T{j} = t_mat(:,1); % only need unique t
    S{j} = sum(s_mat,2); % sum s with same i
    % if T{j} is scalar, then dT{j} is empty
    % dT{j} = diff(T{j}) - mean(diff(T{j})); % fluctuation rel, to mean period
    % increment search
    clear *_mat l
    j = j+1;
end; % while
clear i j k l s t temp1 temp2

% extract composite signal
T = cellextract(T); S = cellextract(S); T=T(:); S=S(:);
[t_sorted,ind] = sort(T); 
s_sorted = S(ind);
num_nan = length(t_sorted) - sum(isnan(t_sorted)); % number of nans
if(num_nan); t_sorted = t_sorted(1:num_nan); s_sorted = s_sorted(1:num_nan); end;
clear ind S T

% determine whether time delay is detectable
% dT  = cellextract(dT); % dT(j,:) = dT{j} = period - mean; jth observation
% ddT = dT(abs(dT) > T_res); % 1 if fluctuation > timing resolution of pulsar

%% Plotting
% colourbar shows impact parameter
thermal_map = thermal2(N_lens); 
parfor j=1:N_lens
    colour_hsv = rgb2hsv(thermal_map(j,:));
    map_hsv = [repmat(colour_hsv(1:2),[N_pulse 1]) linspace(0,colour_hsv(3),N_pulse)'];
    colour{j} = hsv2rgb(map_hsv);
    cbar_map(j,:) = colour{j}(N_pulse,:);
end % for
colour_sorted = colormap_helper([0 0 0;.8 .8 .8],length(t_sorted)); % grey
% get axis handles for subplots
ax = plot_axes((N_lens > 1),2,2,{[1 2],3,4}); 

% plot composite signal
set(ax.figure,'CurrentAxes',ax.handle(1)); hold on;
scatter(t_sorted,s_sorted-1,36,colour_sorted,'Marker','*','DisplayName','signal');
hLine = plot(t_sorted,s_sorted-1,'k-');
set(get(get(hLine,'Annotation'),'LegendInformation'),...
   'IconDisplayStyle','off'); % Exclude line from legend
% plot s(t) for original signal
% scatter(t_source,s_source,36,colour,'Marker','*','DisplayName','source');
% plot s(t) for images: all images from one signal are one colour
% if(N_lens > 1)
%     parfor j=1:N_lens  % different marker
%         hILines(j) = scatter(t_images(:,j),s_images(:,j),36,colour,'Marker',marker{j});
%     end; clear j % for
%     hCGroup = hggroup; set(hILines,'Parent',hCGroup);
%     set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
%     set(hCGroup,'DisplayName','signal from images (per source)');
% end; % if
xlabel('Time (s)','fontsize',10,'fontname','Century Schoolbook L'); 
ylabel('Flux relative to source \mu(t)-1','fontsize',10,'fontname','Century Schoolbook L'); 
title('Signal from pulsar','fontsize',10,'fontname','Century Schoolbook L','fontweight','bold');
axis tight; clear h*
ticklabelformat(ax.handle(1),'xy','%2.6g');

% plot time delay vs emission time
set(ax.figure,'CurrentAxes',ax.handle(2)); hold on;
% cannot display total td(t) since each part of the td is from a different
% image, hence original pulse time is different 
% scatter(t_sorted,ones(size(t_sorted))*1e6,36,colour_sorted,'Marker','*');
parfor j=1:N_lens  % different marker
    scatter(t_source,tau_array(:,j)*1e6,36,colour{j},'Marker','.');
end; clear j % for
xlabel('Emission time (s)','fontsize',10,'fontname','Century Schoolbook L'); 
ylabel('Time delay (\mus)','fontsize',10,'fontname','Century Schoolbook L');
title('Time delay','fontsize',10,'fontname','Century Schoolbook L','fontweight','bold'); 
axis tight; clear h*
ticklabelformat(ax.handle(2),'xy','%2.6g');

%% Output
t_struct.t = t_sorted; t_struct.mu = s_sorted; 
t_struct.ax = ax.handle; t_struct.fig = ax.figure;
t_struct.colour = colour; t_struct.cbar_map = cbar_map;

end % function