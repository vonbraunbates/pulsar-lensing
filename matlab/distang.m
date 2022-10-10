% distang.m
% Calculate angular diameter distances in general Dyer-Roder metric
% INPUT: Omega_0: vector of dimensionless fraction of cosmological fluids
%        z_0: redshift of observer (usually 0)
%        z_1: redshift of source
% OPTIONAL: clumping parameter alpha
%           constant dark energy of state param. w0
%           z-varying parameter w1
% OUTPUT: D = D_ang / (c/H0) dimensionless dist

function DA2 = distang(Omega_0,z_0,z_1,alpha,w0,w1)
hold all;
options = odeset('RelTol',1e-5,'AbsTol',1e-5);
if(sum(Omega_0)~=1)
    disp('Your Omega quantities do not sum to 1! Exiting.');
    return
end

%% Varying inputs for cosmological fluids
switch(nargin) % optional dark energy params
    case(3)
        alpha = ones(size(Omega_0)); w0 = -1; w1 = 0;
    case(4)
        w0 = -1; w1 = 0; 
end

%% Generalised beam eqn.
% Initial conditions:
w_0 = [0 w0+w1*z_0 -1/3];
H_0 = sqrt(sum(Omega_0.*(1+z_0).^(3.*(1+w_0))));
y0 = [(H_0*(1+z_0))^(-1) 0];% [dr/dz r]|z=z_0
[z y] = ode45(@dyer_roder,[z_0 z_1],y0,options); % was tspan = z_0:0.01:z_1
DA2 = y(end,2);
%plot(z,DA2,'c-','DisplayName',['\alpha = ',num2str(alpha(1))]); 

%% Generalised beam eqn. integrator
    function dy = dyer_roder(z,y)
        dy = zeros(2,1); % d/dz [dr/dz r]
        w = [0 w0+w1*z -1/3]; 
        Omega = Omega_0(1:2).*(1+z).^(3.*(1+w(1:2)))./ ...
                (sum(Omega_0.*(1+z).^(3.*(1+w))));
        q = 0.5*sum(Omega_0(1:2).*(1+3.*w(1:2)) ...
            .*(1+z).^(1+3.*w(1:2))) ./ ...
            (sum(Omega_0.*(1+z).^(1+3.*w)));
        dy(1) = -(z+1)^(-1)*((y(1)*(3+q) + ... 
            y(2)/2*(3*sum((1+w(1:2)).*alpha(1:2).*Omega))/(1+z)));
        dy(2) = y(1);
    end

end