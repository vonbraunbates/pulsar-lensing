% Generate observational signatures from the lensing results
function t_struct = plot_tdelay(t_delay,mu,t_lens,T_res)
    
    % get time delays and magnification
    N_lens    = numel(t_delay);         % t is a cell: t{i} = dt(ith lens)
    N_pulse   = size(t_delay{1},1);
    N_images  = max(cellfun(@(x)size(x,2), t_delay));
    tau_array = cellextract(t_delay)';  % t_arr(:,i)  = t{ith lens}(:,:)
    mu_array  = cellextract(mu)';       % mu_arr(:,i) = mu{ith lens}(:,:)
    
    % setup signal
    t_struct.T_res = T_res;               % timing residual (s)
    s_source = ones(size(t_lens));        % original pulsar signal
    
    %% Create signals
    % add lensing effects to signal
    t_images = t_lens + tau_array;
    t_images = reshape(t_images,[N_pulse N_lens*N_images]);
    s_images = s_source .* mu_array;
    s_images = reshape(s_images,[N_pulse N_lens*N_images]);
    clear mu_array
    
    % bin signals by time
    t_images = t_images(:); s_images = s_images(:);
    [temp1,ind] = sort(t_images);  % sort t
    temp2 = s_images(ind);         % sort mu by t value
    temp3 = t_lens(ind);           % sort emission time by t value
    clear ind;
    k = 1; l = 1;
    while k < length(temp1);
        % find signals close together
        dt = temp1 - temp1(k);
        t = temp1((dt >= 0) & (dt < 10*T_res));
        s = temp2((dt >= 0) & (dt < 10*T_res));
        d = temp3((dt >= 0) & (dt < 10*T_res));
        k = find(dt > 10*T_res,1,'first');
        % bin only those signals
        [counts, bin] = histc(t, [min(t) : T_res : max(t)+T_res]);
        max_counts = max(counts); % largest nr of superposed signals
        m=0; % number of non-empty bins
        % Non-empty bins contain signals which will be superimposed
        for i = 1:length(counts)
            if(~isempty(t(bin==i)))
                m = m+1;
                t_mat(m,:) = vec2mat(t(bin==i),max_counts,NaN);
                s_mat(m,:) = vec2mat(s(bin==i),max_counts,0);
                d_mat(m,:) = vec2mat(d(bin==i),max_counts,NaN);
            end % if
        end; % for
        clear bin i max_counts nbins
        % get time, signal for composite
        T{l} = t_mat(:,1);          % only need unique t
        S{l} = sum(s_mat,2);        % sum s with same i
        D{l} = d_mat(:,1);
        % increment search
        clear *_mat dt
        l = l+1;
    end; % while
    
    clear i k l m s t temp1 temp2
    
    % extract composite signal
    T = cellextract(T); T = T(:); [t_sorted,ind] = sort(T);
    S = cellextract(S); S = S(:); s_sorted = S(ind);
    D = cellextract(D); D = D(:); d_sorted = D(ind);
    clear ind D S T
    
    %% Plotting
    % colourbar shows impact parameter
    thermal_map = ...
        [1.0000   0.7857    0.0357
        1.0000    0.5714    0.0714
        0.9857    0.3643    0.1143
        0.9143    0.1857    0.1857
        0.6714    0.0643    0.3714
        0.4000         0    0.5286
        0.1500         0    0.6000];
    thermal_map = colormap_helper(thermal_map, N_lens);
    for j=1:N_lens
        colour_hsv = rgb2hsv(thermal_map(j,:));
        map_hsv = [repmat(colour_hsv(1:2),[N_pulse 1]) linspace(0,colour_hsv(3),N_pulse)'];
        colour{j} = hsv2rgb(map_hsv);
        cbar_map(j,:) = colour{j}(N_pulse,:);
    end % for
    grey_map = repmat([0 .25 .5 .75]',[1 3]);
    colour{j+1} = colormap_helper(grey_map, numel(t_sorted));
    cmap = vertcat(colour{:});              % concatenated maps
    clen = cellfun(@(x)(size(x,1)),colour); % length of each map
    csum = cumsum(clen) - clen;             % starting index of each map
    
    % get axis handles for subplots
    ax = plot_axes((N_lens > 1),2,2,{1,2,3,4}); % {[1 2],3,4});
    
    % Line style default for single lens, dotted for multiple
    if(N_lens==1); linespec = 'k-'; else linespec = 'k:'; end; % if
    
    % plot radio signal from each image
    set(ax.figure,'CurrentAxes',ax.handle(1)); hold on; colormap(cmap);
    ctemp = csum(end) + [1:length(t_sorted)]'; % end was j
    hLine = plot(t_sorted(:), s_sorted(:), linespec);
    set(get(get(hLine,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    scatter(t_sorted(:), s_sorted(:), 360, ctemp, 'Marker','.');
    set(gca,'CLim',[1 sum(clen)]); clear ctemp; freezeColours(gca);
    xlabel('Time (s)','interpreter','latex');
    ylabel('$\mu$','interpreter','latex','rotation',0);
    title('{\bf Signal from pulsar}','interpreter','latex');
    axis tight; clear h*
    ticklabelformat(gca,'xy','%2.6g');
    set(gca,'XTickLabel',get(gca,'Xticklabel'),'FontName','Courier 10 Pitch','FontSize',20,'fontweight','bold');
    set(gca,'XTickMode','auto','XTickLabelMode','auto');
    
    % plot change in time delay vs observation time
    set(ax.figure,'CurrentAxes',ax.handle(2)); colormap(cmap); hold on;
    dt = nan([1 length(t_sorted)]);
    ctemp = csum(end) + [1:length(t_sorted)]';
    dt(2:end) = 1e6*diff(t_sorted - d_sorted, 1, 1);
    scatter(t_sorted, dt, 360, ctemp, 'Marker','.');
    plot(t_sorted, dt,linespec);
    set(gca,'CLim',[1 sum(clen)]); freezeColours(gca);
    xlabel('Emission time (s)','interpreter','latex');
    ylabel('$d\tau$','interpreter','latex','rotation',0);
    title('{\bf Change in time delay}','interpreter','latex');
    axis tight; clear h*
    ticklabelformat(gca,'xy','%2.6g');
    set(gca,'XTickLabel',get(gca,'Xticklabel'),'FontName','Courier 10 Pitch','FontSize',20,'fontweight','bold');
    set(gca,'XTickMode','auto','XTickLabelMode','auto');
    
    % plot time delay per lens
    set(ax.figure,'CurrentAxes',ax.handle(4)); hold on; colormap(cmap);
    for j=1:N_lens
        ctemp = csum(j) + [1:clen(j)]'; % end was j
        if(N_lens~=1) % change in time delay if > 1 lenses
            dt(2:N_pulse,j) = 1e6*diff(tau_array(:,j) - min(tau_array(:,j)), 1, 1);
            scatter(t_lens(:,j), dt(:,j), 360, ctemp, 'Marker','.');
            plot(t_lens(:,j), dt(:,j),linespec);
            ylabel('$d\tau$','interpreter','latex','rotation',0);
            title('{\bf Change in time delay per lens}','interpreter','latex');
        else % time delay if 1 lens
            ttemp = (tau_array - min(tau_array))*1e6;
            scatter(t_lens, ttemp, 360, ctemp, 'Marker','.');
            plot(t_lens, ttemp,linespec);
            ylabel('$\tau$','interpreter','latex','rotation',0);
            title('{\bf Relative time delay per lens}','interpreter','latex');
        end; % if
    end; clear j % for
    set(gca,'CLim',[1 sum(clen)]); freezeColours(gca);
    xlabel('Emission time (s)','interpreter','latex');
    axis tight; clear h*
    ticklabelformat(gca,'xy','%2.6g');
    set(gca,'XTickLabel',get(gca,'Xticklabel'),'FontName','Courier 10 Pitch','FontSize',20,'fontweight','bold');
    set(gca,'XTickMode','auto','XTickLabelMode','auto');
    
    %% Output
    t_struct.t = t_sorted; t_struct.mu = s_sorted;
    t_struct.ax = ax.handle; t_struct.fig = ax.figure;
    t_struct.colour = colour; t_struct.cbar_map = cbar_map;
    
end % function
