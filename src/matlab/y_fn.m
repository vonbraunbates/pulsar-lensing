function y = y_fn(x)
        X = abs(x);
        % 2. mass distribution m(|x|)
        switch(num_model)
            case{1}
                m = ones(size(x));
                
            case{2} % m/x = x/x_max^2; 1/x
                m(X <= x_max) = x(X <= x_max).^2./x_max.^2;
                m(X > x_max)  = 1;
                
            case{3}
                g(X > 1)  = 2./sqrt(X(X > 1).^2 - 1) .* atan(sqrt((X(X > 1) - 1)./(X(X > 1) + 1)));
                g(X < 1)  = 2./sqrt(1 - X(X < 1).^2) .* atanh(sqrt((1 - X(X < 1))./(1 + X(X < 1))));
                g(X == 1) = 1;
                m = 4*kappa_s*(log(X/2) + g);
        end
        % rho_s = f(x) same as 0 = f(x) - rho_s
        y = x - m./x - rho_s;
    end % function