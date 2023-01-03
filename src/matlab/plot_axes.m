
function ax = plot_axes(flag,nrows,ncols,varargin)
    % set figure size, visibility    
    % set(0,'DefaultFigureVisible','off');
    scrsz = get(0,'ScreenSize');
    set(0,'DefaultFigurePosition',[1 .9*scrsz(4) .99*scrsz(3) .9*scrsz(4)]);
        
    % get handle for entire figure
    ax.figure = figure('visible','on');
    n = nrows*ncols;
    
    % axes in plot
    if(flag)
        ax.gap = [0.14 0.08]; % room for colourbar
        ax.max = [1.08 1.0]; 
        ax.min = [0.08 0.08];
    else
        ax.gap = [0.08 0.08]; % no room
        ax.max = [1.04 1.0]; 
        ax.min = [0.08 0.08];
    end
    ax.size = (ax.max - ax.min)./[ncols nrows]; % outer size of axes
    ax.box  = ax.size - ax.gap;                 % inner size of axes
    ax.coord(:,1) = ax.min(1) + ax.size(1).*mod([1:n]-1,ncols); % x starting co-ord
    ax.coord(:,2) = ax.max(2) - ax.size(2).*ceil([1:n]./ncols); % y starting co-ord
    ax.coord(:,3) = ax.box(1);                  % width
    ax.coord(:,4) = ax.box(2);                  % height
    
    % get handles for subplots
    if(nargin == 3)
        for i=1:n
            ax.handle(i) = subplot(nrows,ncols,i,'Parent',ax.figure); 
        end % for
    else
        for i=1:numel(varargin{1}) % subplot cannot be vectorised
            P = sort(varargin{1}{i}(:),'ascend');    
            switch(numel(P)) % test for size of P{i}

                case(1) % P{i} is a scalar
                    
                case(2) % P{i} is a vector
                    if(all(mod(P,ncols) == mod(P(1),ncols))) % P is a column
                        % height is from bottom of P(2) to top of P(1)
                        ax.coord(P(1),4) = sum(ax.coord(P,4)) + ax.gap(2)*floor((P(2)-P(1))/ncols); 
                        % bottom corner is the y co-ord of P(2)
                        ax.coord(P(1),2) = ax.coord(P(2),2);      
                    elseif(all(ceil(P/ncols) == ceil(P(1)/ncols))) % P is a row
                        % width is from left of P(1) to right of P(2)
                        ax.coord(P(1),3) = sum(ax.coord(P([1 2]),3)) + ax.gap(1)*(P(2)-P(1)); 
                    end % if
                    
                case(4) % P is a matrix
                    % width is from left of P(1) to right of P(2)
                    ax.coord(P(1),3) = sum(ax.coord(P([1 2]),3)) + ax.gap(1)*(P(2)-P(1)+1); 
                    % height is from bottom of P(3) to top of P(1)
                    ax.coord(P(1),4) = sum(ax.coord(P([1 3]),4)) + ax.gap(2)*floor((P(3)-P(1))/ncols);
                    % bottom corner is the y co-ord of P(3)
                    ax.coord(P(1),2) = ax.coord(P(3),2); 
                    
            end % switch
            ax.handle(i) = subplot('position',ax.coord(P(1),:),'Parent',ax.figure);
        end % for
    end % if
    
    % Remove warnings in legend
    warning('off','MATLAB:legend:UnsupportedFaceColor');
    warning('off','MATLAB:legend:PlotEmpty');
    warning('off','MATLAB:legend:IgnoringExtraEntries');
    
end % function