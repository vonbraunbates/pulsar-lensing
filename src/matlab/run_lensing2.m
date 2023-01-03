function run_lensing2(filename)
        % read in data
        struct = load(filename); % read in data
        names = fieldnames(struct); % get variable names
        for l=1:numel(names) % assign variable values to names
            eval([names{l},' = getfield(struct,names{l})']);
        end
        clear l struct names
        
        result = lensingdb(dist_cell,M_arr,num_model,scale_mode,b_lens);
        names = fieldnames(result); % get variable names
        for l=1:numel(names) % get variable values from result
            eval([names{l},' = getfield(result,names{l})']);
        end;
        clear l names result
        
        save(filename,'-append');           % save to .mat file
        saveas(handles.fig,filename,'psc2'); % save figure in .eps format
        
        close all; % clear workspace
end % function