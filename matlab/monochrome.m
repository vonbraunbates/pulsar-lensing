% monochrome.m
% Given an rgb colour, create a monochromatic colourmap with n iterations in the colourmap
function colourmap = monochrome(colour,n)
    colour_hsv = rgb2hsv(colour);
%     map_hsv(:,1) = colour_hsv(1)*ones(n,1); % keep hue
%     map_hsv(:,2) = colour_hsv(2)*ones(n,1); % keep sat
%     map_hsv(:,3) = linspace(0.1,1,n); % adjust hue
map_hsv = [repmat(colour_hsv(1:2),[n 1]) linspace(0.1,1,n)'];
map = hsv2rgb(map_hsv);    
colourmap = colormap_helper(map,n);
end % function