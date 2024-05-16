
function [TURB,MOM] = cloud_mask(TURB,MOM,MC,rh_limit,ext_factor,nan_replace)

if nargin<6 || isempty(nan_replace)
    nan_replace = false;
end
if nargin<5
    ext_factor = 0;
end
if nargin<4
    rh_limit = 98;
end


Nseg = size(TURB,1);

for i_s = 1:Nseg
    
    % Interpolate MC mask
    
    mask_interp = interp1(MC(i_s).time,MC(i_s).CLOUD_mask,TURB(i_s).time,'linear');
    TURB(i_s).MC_mask = double(mask_interp>0.5);
    TURB(i_s).MC_mask(isnan(mask_interp)) = nan;
    
    TURB(i_s).MC_mask_ind = mask2ind(TURB(i_s).MC_mask,nan_replace);
    if size(TURB(i_s).MC_mask_ind,1)>0
        TURB(i_s).MC_mask_timestamp = [TURB(i_s).time(TURB(i_s).MC_mask_ind(:,1)),...
            TURB(i_s).time(TURB(i_s).MC_mask_ind(:,2))] + duration(0,0,0.5/TURB(i_s).fsamp);
    else
        TURB(i_s).MC_mask_timestamp = [];
    end 
    
    MOM.MC_nan_fraction(i_s) = sum(isnan(TURB(i_s).MC_mask)) / length(TURB(i_s).MC_mask);
    MOM.MC_cloud_fraction(i_s) = sum(TURB(i_s).MC_mask,'omitnan') / sum(~isnan(TURB(i_s).MC_mask));  
    
    
    % RH mask
    
    TURB(i_s).RH_mask = double(TURB(i_s).RH>rh_limit);      
    TURB(i_s).RH_mask(isnan(TURB(i_s).RH)) = nan;
    
    TURB(i_s).RH_mask_ind = mask2ind(TURB(i_s).RH_mask,nan_replace);
    if size(TURB(i_s).RH_mask_ind,1)>0
        TURB(i_s).RH_mask_timestamp = [TURB(i_s).time(TURB(i_s).RH_mask_ind(:,1)),...
            TURB(i_s).time(TURB(i_s).RH_mask_ind(:,2))] + duration(0,0,0.5/TURB(i_s).fsamp);
    else
        TURB(i_s).RH_mask_timestamp = [];
    end
    
    MOM.RH_nan_fraction(i_s) = sum(isnan(TURB(i_s).RH_mask)) / length(TURB(i_s).RH_mask);
    MOM.RH_cloud_fraction(i_s) = sum(TURB(i_s).RH_mask,'omitnan') / sum(~isnan(TURB(i_s).RH_mask));
    
    
    % OR mask
    
    TURB(i_s).RH_mask_ext = extend_mask(TURB(i_s).RH_mask,ext_factor,'proportional');
    TURB(i_s).MC_mask_ext = extend_mask(TURB(i_s).MC_mask,ext_factor,'proportional');
    
    sum_mask = sum([TURB(i_s).MC_mask_ext,TURB(i_s).RH_mask_ext],2,'omitnan');
    
    TURB(i_s).OR_mask = double(sum_mask>0);
    TURB(i_s).OR_mask(isnan(sum_mask)) = nan;
    
    TURB(i_s).OR_mask_ind = mask2ind(TURB(i_s).OR_mask,nan_replace);    
    if size(TURB(i_s).OR_mask_ind,1)>0
        TURB(i_s).OR_mask_timestamp = [TURB(i_s).time(TURB(i_s).OR_mask_ind(:,1)),...
            TURB(i_s).time(TURB(i_s).OR_mask_ind(:,2))] + duration(0,0,0.5/TURB(i_s).fsamp);
    else
        TURB(i_s).OR_mask_timestamp = [];
    end
    
    MOM.OR_nan_fraction(i_s) = sum(isnan(TURB(i_s).OR_mask)) / length(TURB(i_s).OR_mask);
    MOM.OR_cloud_fraction(i_s) = sum(TURB(i_s).OR_mask,'omitnan') / sum(~isnan(TURB(i_s).OR_mask));
    
end

end