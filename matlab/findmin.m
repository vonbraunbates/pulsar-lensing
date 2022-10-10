function xmin = findmin(name,param_cell,xc,yc) %funstr,arglist)
% INPUT: [xc,yc]: co-ord boundary region in which to find minimum,
%           delimited by NaN
% OUTPUT: xgridmin: location of minima
%         zmin: value of minima
rho_s = param_cell{1}(1); phi_s = param_cell{1}(2);
switch(name)
    case('nfw')
        kappa_s = param_cell{2};
        rmax = param_cell{3};
end

%% Minima-finding function
% Discrete closed contours are delimited by NaN
delim = [1;find(isnan(xc))];
for l=1;%:length(delim)-1
    %flag = 1;
    % treat each contour separately
    % c = [xc(delim(l):delim(l+1)-1), yc(delim(l):delim(l+1)-1)];
    % [cphi crho] = cart2pol(c(:,1),c(:,2)); clear cphi
    % find one point within each contour
    %while flag
    %    xp = rand*(max(c(:,1)) - min(c(:,1))) + min(c(:,1));
    %    yp = rand*(max(c(:,2)) - min(c(:,2))) + min(c(:,2));
    %    if inpolygon(xp,yp,xc,yc)
    %        flag = 0;
    %        [theta r] = cart2pol(xp,yp); clear theta xp yp
            % Use each point to start finding exact minima
            rmin(l) = fzero(@extrematwo,rho_s); 
            % was fminsearch
            
    %    end
    %end
end

[xmin(:,1) xmin(:,2)] = pol2cart(phi_s,rmin);

%% Extrema-finding grid function
    function D = extrematwo(r)
        switch(name)
            case('nfw')
                if(r<1)
                    g = 2./sqrt(1 - r^2)*atanh(sqrt((1 - r)./(1 + r)));
                elseif(r>1)
                    g = 2./sqrt(r^2 - 1).*atan(sqrt((r - 1)./(1 + r)));
                elseif(r==1)
                    g = 1;
                end
                % dimenionless mass m(x) = 4*rhos*rs/sigmacr*g(x)
                kappa_s = kappa_s/(4*pi*(log(1+rmax) - (rmax/(1+rmax))));
                m = 4*kappa_s.*(log(r/2) + g); %*rs*rhos/sigmacr=1;
                
            case('point')
                m = 1; % scaled units
        end % switch
        r_si = r - m./r; % lens eqn. y = x - m(x)/x
        D = r_si - rho_s;
    end % function
end % function