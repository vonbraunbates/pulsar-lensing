T_res = .3;
t_images = []; s_images = []; T = {}; S = {};
parfor k=1:10
    t_images = vertcat(t_images,rand(10,1)+10*k);
end;
s_images = rand(size(t_images));
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
    [counts, bin] = histc(t,[min(t):T_res:max(t)+T_res]);
    max_counts = max(counts); % largest nr of superposed signals
    l=1; % number of bins non-empty bins
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
    dT{j} = diff(T{j}); % period per observation
    ddT{j} = dT{j} - mean(dT{j}); % fluctuation rel, to mean period
    % increment search
    clear *_mat l
    j = j+1;
end; % while
clear i j k l s t temp1 temp2

T = cellextract(T); S = cellextract(S);
[t_sorted,ind] = sort(T(:));
s_sorted = S(:); s_sorted = s_sorted(ind);