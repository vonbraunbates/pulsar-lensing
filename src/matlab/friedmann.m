function [z DM] = friedmann(zspan,Omega,w0,w1,gamma,h)

% Initial conditions into single matrix
Omega_m = Omega;
Omega_x = 1 - Omega_m;
y0 = [Omega_m Omega_x 0.0];

% Solve coupled ODEs, integrating from z=0
options = odeset('RelTol',1e-5,'AbsTol',1e-5);
[z,y] = ode45(@friedmann_rhs,[0;zspan],y0,options);

% Calculate DM = 5log[(c/H0/10pc)*Int_0^z{H0/H(z')dz'}]
E_int = y(:,3);
% 43.089 = -5log10(c/H0/10 pc)
DM = 5*log10((1+z).*E_int) + 43.0989 - 5*log10(h);

    function dy = friedmann_rhs(z,y)
        dy = zeros(3,1);    % a column vector
        H = sqrt(y(1)+y(2)); % H(z)
        dy(1) = (1+z)^(-1) * (3*y(1) - gamma/H); %dOmega_m/dz
        dy(2) = (1+z)^(-1) * (3*y(2)*(1+(w0 + w1*z/(1+z))) + gamma/H); %dOmega_x/dz
        dy(3) = 1/H; % Integrand in luminosity distance
    end
end
