% lens_contour.m
% Create contour plot of the arrival time of null geodesics from source to
% observer via an elliptical gravitational lens
close all
%% Insert parameters - can be vectorised
ask = input('Do you want to change current parameters? y/n/d (default): ','s');
if(ask=='y')
    %    x = input('Input position on lens plane: [x1 x2]');
    b = input('Deflection scale of potential in (0 1): ');
    s = input('Dispersion of mass from centre in [0 1]');
    e = input('Ellipticity of lens galaxy in [0 1]: ');
    y = input('Input position on source plane: [y1 y2]');
elseif(ask=='n')
    disp('Fine.');
elseif(ask=='d')
    b = [0:0.25:1]; % b = 4*pi*D_LS/D_OS*(sigma/c)^2 << 1 ??
    s = [0:0.25:1]; % s probably closer to zero
    e=0.1;%e = [0:0.2:1]; % eccentricity
    y = [1 1];
    disp('Using arrays.');
else
    disp('Error: please enter y/n/d!');
end
% Plane of the lens: lens at origin
[x1,x2] = meshgrid(-10:.5:10, -10:.5:10);

%% Calculate time delay
% Geometric time delay
tau_geom = 0.5*((x1-y(1)).^2+(x2-y(2)).^2);
for k=1:length(e)
    %figure(k)
    for j=1:length(s)
        figure(j)
        for i=1:length(b)
            % Gravitational time delay
            tau_pot = b(i)*sqrt(s(j)^2 + (1-e(k))*x1.^2 + (1+e(k))*x2.^2);
            % Total
            tau = tau_geom - tau_pot;
            
            %% Plot the corresponding maps
            subplot(ceil(length(b)/3),1+ceil(length(b)/3),i);
            
            meshc(x1, x2, tau_geom);%,'FaceColor','none');
            colormap('gray'); freezeColours; hold on;
            surfc(x1,x2,tau); colormap('jet'); shading interp;
            hold off;
            title(['b = ',num2str(b(i))]); %,' s = ',num2str(s(j))]);
            set(gca,'ztick',[])
            axis square;
        end % for
        %end
        set(gca,'xtickMode', 'auto')
        % title for entire figure
        mtit(['Lens plane contours: e = ',num2str(e(k)),' s = ',num2str(s(j)),' source at ',mat2str(y)], ...
            'FontWeight','bold','FontSize',14);
        % create colourbar for each figure
        h=findobj(gcf,'type','axes');
        q=get(h,'clim');
        set(h,'clim',[min([q{:}]) max([q{:}])]);
        axes('Position', [0.05 0.05 0.95 0.9], 'Visible', 'off');
        set(gca,'clim',[min([q{:}]) max([q{:}])]); % MUST go before
        cbar = colorbar('location','EastOutside'); % THIS line
        % export to PNG file
        figname = ['lens_contour_s',num2str(10*s(j)),'_src',mat2str(y)];
        %figname = ['lens_contour_e',num2str(10*e(k)),'_src',mat2str(y)];
        print('-dpng',figname);
    end % for
end % for