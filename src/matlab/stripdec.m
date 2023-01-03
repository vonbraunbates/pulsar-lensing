function str_cell = stripdec(num_array,dec_char)

% if num is a matrix, separate entries
num_cell = num2cell(num_array);
parfor i=1:numel(num_cell)
    num = num_cell{i};
    % convert number to string
    sprintf_num = sprintf('%0e',num);
    % find location of decimal point
    dec_regexp = strcat('\',dec_char);
    dec_ind = regexpi(sprintf_num,dec_regexp);
    % find exponent value
    exp_ind = regexpi(sprintf_num,'[\+\-]')-1; % location of e,E
    exp_val = str2num(sprintf_num(exp_ind+1:end));
    
    % split number before, after decimal point
    sprintf_ante = sprintf_num(1:dec_ind-1);
    sprintf_post = sprintf_num(dec_ind+1:exp_ind-1);
    % find index of last preceeding zero
    preceed_ind = regexpi(sprintf_ante,'(?<=[1-9]0+\.)','end');
    if(isempty(preceed_ind)); preceed_ind=0; end;
    % find index of first trailing zero
    trail_ind = regexpi(sprintf_post,'0+$','start');
    if(isempty(trail_ind)); trail_ind=length(sprintf_post)+1; end;
    % remove preceeding, trailing zeros
    sprintf_ante = sprintf_ante(preceed_ind+1:end);
    sprintf_post = sprintf_post(1:trail_ind-1);
    % shift to move decimal point after non-zero nos
    exp_shift = dec_ind - trail_ind-1;
    % adjust exponent to retain same number
    exp_val = exp_val + exp_shift;
    if(exp_val < 0)
        exp_sign = '-';
        exp_val = abs(exp_val);
    else
        exp_sign = '+';
    end
    
    % output number as string
    str_mantissa = strcat(sprintf_ante,sprintf_post);
    str_exp = strcat(exp_sign,num2str(exp_val));
    str = strcat(str_mantissa,'e',str_exp);
    
    str_cell{i} = str;
    
end % for
end % function