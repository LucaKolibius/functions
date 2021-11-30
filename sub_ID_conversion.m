function [ conv_sub_ID ] = sub_ID_conversion( sub_ID, has_micro )

% if( contains(sub_ID, 'P') == 1 || contains(sub_ID, 'ERL') == 1 ) % does not work on R2015b
    
if strcmp(sub_ID(1), 'P') || any(regexp(sub_ID, 'ERL'))
    %% old to new
    if(contains(sub_ID, 'ERL') == 1)
        x1 = '1';
    elseif(contains(sub_ID, 'P') == 1 )
        x1 = '0';
    end
    if has_micro == 1
        x2 = '0';
    else; x2 = '1';
    end
    conv_sub_ID = ['sub-' x1 x2 sub_ID(2:3)];
% elseif( contains(sub_ID, 'sub') == 1 ) % see above
elseif any(regexp(sub_ID, 'sub'))
    %% new to old
    if( str2double(sub_ID(5)) == 1 )
        x1 = 'ERL';
    elseif( str2double(sub_ID(5)) == 0 )
        x1 = '';
    end
    conv_sub_ID = ['P' sub_ID(end-1:end) x1];
end


end

