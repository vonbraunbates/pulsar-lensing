function x_0 = rootsearch(f,df,d2f,a,b)

eps_ = double(eps('single'));
d2x0 = []; dx0 = []; x0 = []; x_0 = [];
options = optimset('FunValCheck','on', ... % f(x0) finite
    'TolFun',eps_/1e1); % tolerance f(x)

%%{
% Plot functions to check:
xx = linspace(a,b,1e4);
figure(6); hold all
plot(xx,zeros(1e4,1),'k-','DisplayName','f = 0','MarkerSize',4)
plot(xx,d2f(xx),'.','DisplayName','d^2f/dx^2','MarkerSize',4)
plot(xx,df(xx),'.','DisplayName','df/dx','MarkerSize',4)
plot(xx,f(xx),'.','DisplayName','f(x)','MarkerSize',4)
axis([a b -10 10]); legend('show','Location','Best');
%}
    
% find where d2f/dx2 changes sign:
if (sign(d2f(a)) ~= sign(d2f(b)))
    d2x0 = fzero(d2f,[a b],options);
else
    d2x0 = NaN;
end

d2x0 = d2x0(isfinite(d2x0)); % remove NaN
d2x0(abs(d2x0) < eps_) = 0;

% between each root look for roots of df/dx
rangeint = unique([a;d2x0;b]); % sort

for i = 2:length(rangeint)
    if(~isfinite(df(rangeint(i-1)))) % df(a) = +/- Inf
        if (sign(df(rangeint(i-1) - eps_)) ~= sign(df(rangeint(i))))
            int = [rangeint(i-1) - eps_ rangeint(i)];
            dx0(i-1) = fzero(df,int,options);
        elseif (sign(df(rangeint(i-1) + eps_)) ~= sign(df(rangeint(i))))
            int = [rangeint(i-1) + eps_ rangeint(i)];
            dx0(i-1) = fzero(df,int,options);
        end
    elseif(~isfinite(df(rangeint(i)))) % df(b) = +/-Inf
        if (sign(df(rangeint(i-1))) ~= sign(df(rangeint(i) - eps_)))
            int = [rangeint(i-1) rangeint(i) - eps_];
            dx0(i-1) = fzero(df,int,options);
        elseif (sign(df(rangeint(i-1))) ~= sign(df(rangeint(i) + eps_)))
            int = [rangeint(i-1) rangeint(i) + eps_];
            dx0(i-1) = fzero(df,int,options);
        end
    elseif(sign(df(rangeint(i-1))) ~= sign(df(rangeint(i))))
        % f continuous over int
        int = rangeint(i-1:i);
        dx0(i-1) = fzero(df,int,options);
    else % df is undefined at a or b
        dx0(i-1) = NaN;
    end
end

dx0(abs(dx0) < eps_) = 0;
dx0 = dx0(isfinite(dx0)); % remove NaN

% between those roots look for roots of f
rangeint = unique([a;d2x0';dx0';b]); % sort

for i = 2:length(rangeint)
    % f = +/- Inf at a or b breaks fzero(f,[a b])
    if(~isfinite(f(rangeint(i-1)))) % f(a) = +/- Inf
        if (sign(f(rangeint(i-1) - eps_)) ~= sign(f(rangeint(i)))) 
            int = [rangeint(i-1) - eps_ rangeint(i)];
        elseif (sign(f(rangeint(i-1) + eps_)) ~= sign(f(rangeint(i)))) 
            int = [rangeint(i-1) + eps_ rangeint(i)];
        end
    elseif(~isfinite(f(rangeint(i)))) % f(b) = +/- Inf
        if (sign(f(rangeint(i-1))) ~= sign(f(rangeint(i) - eps_))) 
            int = [rangeint(i-1) rangeint(i) - eps_];
        elseif (sign(f(rangeint(i-1))) ~= sign(f(rangeint(i) + eps_))) 
            int = [rangeint(i-1) rangeint(i) + eps_];
        end
    elseif(sign(f(rangeint(i-1))) ~= sign(f(rangeint(i)))) 
        % f continuous over int
        int = rangeint(i-1:i);
    end
    % Now find root within modified interval
    try
        [x0(i-1),~,exitflag,~] = fzero(f,int,options);
        if(exitflag==1) 
        x_0 = [x_0;x0(i-1)]; 
        else
        % f is undefined over [a,b]
        x0(i-1) = NaN;
        end
    catch
        % disp('FZERO error.')
    end; % try
end

% concatenate zeros
xvals = [d2x0 dx0];
fvals = f(xvals); % check roots of f', f" zeros of f
x_0 = sort([x_0; xvals(abs(fvals) < eps_)']); % keep true zeros
if(isempty(x_0)); 
    disp('No roots!'); x_0 = NaN;
else
    x_0 = x_0(logical([1,(diff(x_0) > eps_)'])); % remove elements equal within tol
end; % if

% The value x returned by fzero is near 
% a point where fun changes sign, 
% or NaN if the search fails.
% ONLY a zero if fun is continuous
% Otherwise a divergent discontinuity.

plot(x_0,zeros(size(x_0)),'x','DisplayName','zeros')
hold off;

end % function