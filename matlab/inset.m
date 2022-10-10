function [h_main, h_inset]=inset(main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% An example:
% fig1=figure;
% plot(vv,ht(log(rr)));
% fig2=figure;
% plot(vv,ht(log(rr))); axis([0 max(vv) -1e4 1e4]);
% [h_main, h_inset]=inset(fig1, fig2);
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.7;
figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
% change the position of the inset on the main plot
set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size inset_size])
% legend can be added afterwards to the inset
% legend(h_main,'toggle') or legend(h_inset,'toggle')