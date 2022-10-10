function H = hankel_matrix(ord, R, N, varargin)
%{
HANKEL_MATRIX: Generates data to use for Hankel Transforms

	The algorithm used is that from:
		"Computation of quasi-discrete Hankel transforms of the integer
		order for propagating optical wave fields"
		Manuel Guizar-Sicairos and Julio C. Guitierrez-Vega
		J. Opt. Soc. Am. A 21(1) 53-58 (2004)

paper defn: (eqn 1)
H[f(r)]   \equiv 2*pi \int dr f(r)J_p(2*pi*kr)r 
H-1[F(k)] \equiv 2*pi \int dk F(k)J_p(2*pi*kr)k 
RHB defn.:
H[f(r)]   \equiv \int dr f(r)J_p(kr)r 
H-1[F(k)] \equiv \int dk F(k)J_p(kr)k 
scaling:
forward:   H[f] = (T * (f.*s_HT.JR) ) ./ s_HT.JV   (eqn 6a)
backward: ~H[F] = (T * (F.*s_HT.JV) ) ./ s_HT.JR   (eqn 6b)
%}
if(~isempty(varargin))
    fhandle1 = figure('visible','on');
    %fhandle2 = figure('visible','on');
    % set figure size, visibility
    set(0,'DefaultFigureVisible','off');
    scrsz = get(0,'ScreenSize');
    set(0,'DefaultFigurePosition',[1 .9*scrsz(4) .99*scrsz(3) .9*scrsz(4)]);
    % axes in plot
    ncols = 2;     nrows = 2;     len=nrows*ncols;
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
    
    clear len scrsz ax.min ax.max ax.size ax.box
end % if

%% Transformation matrix
if(~isinteger(N)); N = floor(N); end; % int nr of Bessel roots
%	Calculate N+1 roots:
c = bessel_zeros('J',ord,N+1);
% [jn,jm] = meshgrid(c(1:N),c(1:N)); % alpha_{p,1:N}
% Jn = besselj(ord+1,jn); Jm = Jn';
% But meshgrid runs out of memory!
J = besselj(ord+1,c(1:N)'); 
Jn = abs(repmat(J,N,1)); % rows of Jn are copies of J
% Calculate hankel matrix
% [b,~,n] = unique(c(1:N)*c(1:N)'/c(N+1));
% B = besselj(ord,b); clear b
% bessel = reshape(B(n),[N N]); clear B n
C = (2/c(N+1))*besselj(ord,(c(1:N)*c(1:N)')/c(N+1))./(Jn.*Jn'); %c*c' = jn.*jm
clear Jn 

% Co-ordinate vectors: f_n = f(j_n/V); F_m = F(j_m/R); 
V = c(N+1)/R;           % Maximum frequency
r = c(1:N)/V;           % /V instead of *R/c(N+1);   % Radius vector
v = c(1:N)/R;           % Frequency vector

% Scaling: f_qdht = f(x)/m1; F_qdht = F(k)/m2
% F(k) =  ht[f_qdht] * m2 = (C * (f(x)/m1)) * m2;
% f(x) = iht[F_qdht] * m1 = (C * (F(k)/m2)) * m1;
m1 = abs(J')/R;         %% m1 prepares input vector for transformation
m2 = abs(J')/V;         %% m2 prepares output vector for display

%% Analytical soln if necessary
if(~isempty(varargin))
    % input
    f = [];
    
    % transform and inverse transform
    ht  = @(f) (C*(f(:)./m1)).*m2;
    iht = @(F) (C*(F(:)./m2)).*m1;
    for j=1:2
        f2(:,j)   = ht ( f(:,j));  % forward
        fiht(:,j) = iht(f2(:,j));  % backward
    end % for
    clear j
    f2(:,3)   = 2*pi*f2(:,1).*f2(:,2); % convolution thm.
    fiht(:,3) = iht(f2(:,3));
    
    %% Plotting
    title_str = sprintf('N = %8.0g, R_{max} = %8.0g, N/R = %8.0g',[N,R,N/R])
    
    % actual plots
    figure(fhandle1), subplot(1,2,1), hold all,
    plot(r,fiht,'o'), axis tight;
    xlabel('r'), ylabel('f(r)'); 
    
    subplot(1,2,2), hold all, 
    plot(v,f2,'o'), axis tight;
    xlabel('v'), ylabel('F(v)'); 
    
    mtit(title_str);
end % if

%% assign to struct
H = struct('C',C,'r',r,'v',v,'m1',m1,'m2',m2);
clear C r v m1 m2 f f2 fiht

end % function