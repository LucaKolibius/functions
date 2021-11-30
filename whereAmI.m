function prePath = whereAmI()
% isHome = 0; % using Glasgow Local Machine
% isHome = 1; % using my own labtop
% isHome = 2; % using Glasgow Virtual Machine
global prePath

if isdir('//analyse4.psy.gla.ac.uk/project0309/')
    prePath = '//analyse4.psy.gla.ac.uk/project0309/';
elseif isdir('X:/')
    prePath = 'X:/';
elseif isdir('/analyse/Project0309/')
    prePath = '/analyse/Project0309/';
end
    
% switch isHome
%     case 0
%         prePath = '//analyse4.psy.gla.ac.uk/project0309/';
%     case 1
%         prePath = 'X:/';
%     case 2
%         prePath = '/analyse/Project0309/';
% end

end
