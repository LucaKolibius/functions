function whereAmI(isHome)
% isHome = 1; % using my own labtop
% isHome = 0; % using Glasgow Local Machine
global prePath
switch isHome
    case 0
        prePath = '\\analyse4.psy.gla.ac.uk\project0309\';
    case 1
        prePath = 'X:\';
end

end
