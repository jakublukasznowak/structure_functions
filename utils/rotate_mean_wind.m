function MOM = rotate_mean_wind(MOM)
    
ang = 90 - MOM.MEAN_WDIR;
MOM.U = - MOM.MEAN_WSPD .* cosd(ang);
MOM.V = - MOM.MEAN_WSPD .* sind(ang);

ang = 90-MOM.MEAN_THDG;

MOM.UX =  cosd(ang) .* MOM.U + sind(ang) .* MOM.V;
MOM.VY = -sind(ang) .* MOM.U + cosd(ang) .* MOM.V;

end