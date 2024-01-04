
function [TURB,MOM] = cloud_mask(TURB,MOM,PMA,rh_limit,ext_factor)

if nargin<5
    ext_factor = 1;
end
if nargin<4
    rh_limit = 98;
end


Nseg = size(TURB,1);

for i_s = 1:Nseg
    
    % Interpolate PMA mask
    
    mask_interp = interp1(PMA(i_s).time,PMA(i_s).CLOUD_mask,TURB(i_s).time,'linear');
    TURB(i_s).PMA_mask = double(mask_interp>0.5);
    TURB(i_s).PMA_mask(isnan(mask_interp)) = nan;
    
    MOM.PMA_mask_nan_fraction(i_s) = sum(isnan(TURB(i_s).PMA_mask)) / length(TURB(i_s).PMA_mask);
    MOM.PMA_mask_cloud_fraction(i_s) = sum(TURB(i_s).PMA_mask,'omitnan') / sum(~isnan(TURB(i_s).PMA_mask));  
    
    
    % RH mask
    
    TURB(i_s).RH_mask = double(TURB(i_s).RH>rh_limit);      
    TURB(i_s).RH_mask(isnan(TURB(i_s).RH)) = nan;
    
    MOM.RH_mask_nan_fraction(i_s) = sum(isnan(TURB(i_s).RH_mask)) / length(TURB(i_s).RH_mask);
    MOM.RH_mask_cloud_fraction(i_s) = sum(TURB(i_s).RH_mask,'omitnan') / sum(~isnan(TURB(i_s).RH_mask));
    
    
    % OR mask
    
    TURB(i_s).RH_mask_ext  = extend_mask(TURB(i_s).RH_mask, ext_factor,'proportional');
    TURB(i_s).PMA_mask_ext = extend_mask(TURB(i_s).PMA_mask,ext_factor,'proportional');
    
    sum_mask = sum([TURB(i_s).PMA_mask_ext,TURB(i_s).RH_mask_ext],2,'omitnan');
    
    TURB(i_s).OR_mask_ext = (sum_mask>0);
    TURB(i_s).OR_mask_ext(isnan(sum_mask)) = true; % Assume cloud if both masks are nan
    
end

end