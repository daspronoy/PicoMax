function [label] = label(o,option)

if nargin<2
    option = 0;
end

if option
    label = cellfun(@(x)strrep(x,'G','$\Gamma$'),regexp(o.num.QPATH,'([^ ,:]*)','tokens'));
else
    label = cellfun(@(x)strrep(x,'G','\Gamma'),regexp(o.num.QPATH,'([^ ,:]*)','tokens'));
end

