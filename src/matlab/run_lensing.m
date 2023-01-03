% run_lensing.m
% INPUT: filename = file with parameters in it
% OUTPUT: appends lens data to file "filename"
%         creates graphics files "graphics.fig"

% variables
N_lens = [81 8 1]; % nr lenses
num_model = 1; scale_mode = 0;
b_max = [1 3 10 30 100]; % max impact parameter
M_lens = [1e4 1e5 1e6]; % mass per lens
var_dist = input('Vary distance? y/n: ','s');
if(var_dist=='y')
    r = linspace(1e-2,1-1e-2,1e2);
    pdf_r = (r.^2); pdf_r = pdf_r./max(pdf_r); % 1/r^2
    cdf_r = cumsum(pdf_r)./sum(pdf_r); % normalise cdf
    u = rand([1 N_lens(1)]); % uniform dist
    parfor k=1:numel(u)
        [~,ind(k)] = min(abs(cdf_r - u(k))); % find closest cdf_r
    end; clear k
end % if Dd

% loop over varying params
for i=1%:numel(b_max)
    b = b_max(i)*rand([1 N_lens(1)]); % select random b
    if(var_dist=='y')
        Dd = 1e4*r(ind); % select random Dd
    elseif(var_dist=='n')
        Dd = 5e3*ones([1 N_lens(1)]); % have fixed Dd
    end % if Dd
    for j=1:numel(N_lens)
        M_arr = M_lens(j)*ones([1 N_lens(j)]); % M
        if(j==1); b_lens = b; dist_cell = {'dist',sort(Dd),1e4};
        else
            q = randperm(N_lens(1)); % rand permuation of indices
            b_lens = b(q(1:N_lens(j))); dist_cell = {'dist',sort(Dd(q(1:N_lens(j)))),1e4};
        end % if
        b_lens(end+1) = b_max(i); % store b_max
        filename = ['b_',num2str(b_max(i)),'_N_',num2str(N_lens(j)),'_Tres_1e-6'];
        save(filename,'M_arr','dist_cell','b_lens','num_model','scale_mode');
        run_lensing2(filename); % run code
    end % for
end % for
