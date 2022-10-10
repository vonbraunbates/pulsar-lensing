% werner_scaling.m
% Check that results from multiple Einstein rings agree with paper.
% INPUT: M = [M1 M2] in Earth masses
%        dist_mode = {'dist',D1,D2,Ds}

function werner_scaling(M,dist_mode)

% mks units
c = 3E8; G = 6.67E-11; pc = 3.086E16;
M_sol = 1.99e30; M_earth = 5.972e24;
eps_ = eps('single');
M_scale = 4*G*(M_earth/c^2); % or M_sol if M given in solar masses

% Single lens case:
if(nnz(M) == 2) % distances to L(1), L(2)
    D1 = dist_mode{2}*pc; D2 = dist_mode{3}*pc; 
elseif(M(1) == 0) % only S, L(2): set L(1) = L(2) 
    D1 = dist_mode{3}*pc; D2 = dist_mode{3}*pc; 
elseif(M(2) == 0) % only S, L(1): set L(2) = L(1)
    D1 = dist_mode{2}*pc; D2 = dist_mode{2}*pc;
end % if 
Ds = dist_mode{4}*pc; % source
D1s = Ds - D1; D2s = Ds - D2; D12 = D2 - D1; % dist between pairs
% scaling in paper
mu(1)  = M(1)*M_scale*D1s/D1; mu(2) = M(2)*M_scale*D2s/D2;
beta   = (D12*Ds)/(D1s*D2); % relative distance between lenses
zeta_0 = sqrt(sum(mu))*sqrt(D1^2/Ds);
m = mu./sum(mu); % relative lens masses 0 < m < 1

% equations in paper in POLAR co-cordinates (x,phi)
x_E(1) = 0.5*(1 + beta*m(1) + sqrt((1 + beta*m(1))^2 - 4*beta*(m(1)^2)));
x_E(2) = 0.5*(1 + beta*m(1) - sqrt((1 + beta*m(1))^2 - 4*beta*(m(1)^2)));
if((beta == 1 && m(1) == 1) || (beta == 0 || m(1) == 0))
    disp('There should be 1 Einstein ring')
else
    disp('There should be 2 Einstein rings')
end;

% own results
result = multiple_lenses(dist_mode,M',1,[0 0;0 0],[0 0]);
x1 = result.loc_im{:}; rho = hypot(x1(:,1),x1(:,2)); % L(1) images

% scaling
if(m(2)==1) % images are in L(2)
    rho_0 = sqrt(M(2)*M_scale * D2*D2s/Ds);
else % images are in L(1)
    rho_0 = sqrt(M(1)*M_scale * D1*D1s/Ds);
end
r_E = sqrt(x_E)*zeta_0; % eqn for r_E gives the squares of radii
rho_im = unique(rho).*rho_0;

% comparison
disp(['Werner radii are: ',mat2str(r_E,8),' m']);
disp(['Own radii are: ',mat2str(rho_im,8),' m']);

end % function