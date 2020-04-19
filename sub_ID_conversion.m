function [ conv_sub_ID ] = sub_ID_conversion( sub_ID, has_micro )

if( contains(sub_ID, 'P') == 1 || contains(sub_ID, 'ERL') == 1 )
    %% old to new
    if(contains(sub_ID, 'ERL') == 1)
        x1 = '1';
    elseif(contains(sub_ID, 'P') == 1 )
        x1 = '0';
    end
    if(strcmp(has_micro, 'yes') == 1)
        x2 = '0';
    else; x2 = '1';
    end
    conv_sub_ID = ['sub-' x1 x2 sub_ID(2:3)];
elseif( contains(sub_ID, 'sub') == 1 )
    %% new to old
    if( str2double(sub_ID(5)) == 1 )
        x1 = 'ERL';
    elseif( str2double(sub_ID(5)) == 0 )
        x1 = '';
    end
    conv_sub_ID = ['P' sub_ID(end-1:end) x1];
end


end

