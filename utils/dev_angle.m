
function [dir2,dir4] = dev_angle (thdg,wdir)

dirs = ["up","left","down","right","up"];

dev = mod(thdg-wdir,360);

L = length(thdg);
dir4 = repmat("",L,1);
dir2 = repmat("",L,1);

for i = 1:L

    [~,ind] = min(abs( dev(i) - [0 90 180 270 360] ));

    dir4(i) = dirs(ind);

    if ismember(dir4(i),["up","down"])
        dir2(i) = "along";
    elseif ismember(dir4(i),["left","right"])
        dir2(i) = "cross";
    end

end

end