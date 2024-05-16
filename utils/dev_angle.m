
% Determine the orientation of flight segment with respect to mean wind
% using the deviation angle which is the difference between true heading
% and wind direction.
% DIR4 can be up/left/down/right, depending on whichever of 0/90/180/270 is
% closest to the deviation angle.
% DIR2 can be along/across, depending on DIR4 (up/down is along, left/right
% is across)


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