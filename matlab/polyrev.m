function polyrev = polyrev()

%clear all; close all; clc;
grid_res = 0.05; revolve_res=0.2;
x=0:grid_res:5*pi;
z=sin(x+(pi/2));
[polyrev.X,polyrev.Y,polyrev.Z]=polyrevolve(x,z,revolve_res);

%minX = min(min(polyrev.X)); maxX = max(max(polyrev.X)); 
%minY = min(min(polyrev.Y)); maxY = max(max(polyrev.Y));
%[polyrev.XI,polyrev.YI] = meshgrid(minX:grid_res:maxX,minY:grid_res:maxY);
[polyrev.XI,polyrev.YI] = meshgrid(polyrev.X,polyrev.Y);
polyrev.F = TriScatteredInterp(polyrev.X,polyrev.Y,polyrev.Z); 
polyrev.ZI = polyrev.F(polyrev.XI,polyrev.YI);

function [X,Y,Z]=polyrevolve(x,z,n)

% function [X,Y,Z]=polyrevolve(x,z,n)
% ------------------------------------------------------------------------
% This function revolves a 2D polygon around the Z-axis.
%
% x is a vector with x coordinates
% z is a vector with z coordinates
% n is the resolution of the rotation in radians
%
% The vectors x and z are converted to spherical coordinates (rho,r). The
% points are then 'copied' around the z-axis. The space between each point
% is defined by n.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 14/08/2008
% ------------------------------------------------------------------------

[theta, rho, r]=cart2sph(x,zeros(size(x)),z);
THETA=[];
RHO=[];
R=[];
parfor i= 1: length(r)
    if n>(2*x(i))
        %The distance between points on a circle can not exceed the
        %diameter. These points are not revolved.
        %disp(['Warning: n>(2*r(i)), n exceeds maximum value for point no ',num2str(i) ,', this point will not be revolved! '])
        if x(i)==0 %If distance to rotational axis keep point
            theta_step=0;
            THETA=[THETA; theta_step'];
            RHO=[RHO; (rho(i)*ones(1,length(theta_step)))'];
            R=[R; (r(i)*ones(1,length(theta_step)))'];
        end
    else
        theta_inc = real(2*asin((n/2)/x(i))); %Angular spacing between points in radians
        no_steps=round((2*pi)/theta_inc)+1; 
        theta_step=linspace(0,(2*pi),no_steps); 
        theta_step=theta_step(1:end-1);
        THETA=[THETA; theta_step'];
        RHO=[RHO; (rho(i)*ones(1,length(theta_step)))'];
        R=[R; (r(i)*ones(1,length(theta_step)))'];
    end
end

[X Y Z]=sph2cart(THETA,RHO,R);
end % function

end % function