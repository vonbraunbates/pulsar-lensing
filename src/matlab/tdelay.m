function varargout = tdelay(mode,varargin)
% mode{1}: 'z' redshift, 'dist' distance
% mode{2:end}: redshifts or distances, {2...k}:lens, {end}:source
% varargin: arguments if converting time delay

global c t_H pc % mks universal constants

 % dimensionless angular diameter distances
 z = zeros(size(mode)); D = zeros(length(mode));
 switch(mode{1})
    case('z')       % mode redshifts; 
        z(1) = 0; z(2:end) = cell2mat(mode(2:end));
        for i=1:length(mode)
            for j=2:length(mode);
                if( j > i)
                    D(i,j) = distang([0.27 0.73 0],z(i),z(j));
                end % while
            end; clear j % for
        end; clear i % for
        % convert from units of D_H to m
        D = D*c*t_H;

    case('dist')    % mode dist in pc; 
        z(1) = 0; z(2:end) = cell2mat(mode(2:end)).*pc./(c*t_H); 
        % S = C(a:b) returns subcell S of C; 
        % S = C{a,b} returns [C{a}] ... [C{b}]
        D(1,:) = [0 cell2mat(mode(2:end))];
        for i=2:length(mode)
            for j=2:length(mode);
                if( j > i)
                    D(i,j) = abs(mode{i} - mode{j});
                end % while
            end; clear j % for
        end; clear i % for
        % convert from pc to m
        D = D*pc; 
 end % switch
 
 % If not scaling, export distances
 if(nargin == 1)
     D = D(1:size(D,1)-1,2:size(D,2)); z = z(2:end);
     varargout = [{z D}]; 
 % If scaling, unpack varargin, then scale time delays:
 elseif(nargin > 1)
     taumin = varargin{1}; % time delays
     zeta_0 = varargin{2}; % scaling in lens plane
     % Use Dd, Dds, Ds
     Dd = D(1,2); Ds = D(1,3); Dds = D(2,3); zd = z(2);
     % T(x,y) eq. 5.45 Schneider for td at same y (s)
     t_rel = zeta_0^2/c * Ds/(Dd*Dds) * (1+zd).*taumin;
     varargout = {t_rel};
 end % if
end % function