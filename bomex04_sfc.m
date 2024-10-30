%% Settings

% Compute uncertainties?
% (takes quite a long time)
if_uncertainties = false;


% Cloud mask
mask_ql_threshold = 0;    % Relative humidity threshold for cloud mask
mask_ext_factor   = 1;    % Mask extension factor to account for cloud shells


% Structure functions
r_max = 1000; % [m] Max r separation


% Dissipation rate
edr_fit_range = [50 150]; % [m]



%% Load

addpath(genpath(myprojectpath))

load([myprojectpath,filesep,'bomex_turb.mat'])

plotpath = [myprojectpath,filesep,'figures_bomex',filesep,'results'];
if ~isfolder(plotpath), mkdir(plotpath), end


dr = dx; % [m]
seg_length = dx*length(x); % [m]
r_maxlag = r_max/dr;


%% Wind rotation

disp('Wind rotation ...')

Nseg = size(TURB,1);
for i_s = 1:Nseg
    
    TURB(i_s).alt = z(TURB(i_s).ind_z);
    TURB(i_s).dr = dr;
    
    if ~isnan(TURB(i_s).ind_x)
        TURB(i_s).ux =  TURB(i_s).v;
        TURB(i_s).vy = -TURB(i_s).u;
    elseif ~isnan(TURB(i_s).ind_y)
        TURB(i_s).ux =  TURB(i_s).u;
        TURB(i_s).vy =  TURB(i_s).v;
    end
    
    TURB(i_s).ux = detrend(TURB(i_s).ux);
    TURB(i_s).vy = detrend(TURB(i_s).vy);
    TURB(i_s).w = detrend(TURB(i_s).w);
    
    TURB(i_s).pi = TURB(i_s).p/TURB(i_s).rho;
    
end

    
    
%% Thermodynamics

disp('Thermodynamics ...')

% Constants
Rd = 287.04;  % [J/kg/K]
Rv = 461.5;
ep = Rd/Rv;
g  = 9.81;       % [m/s2]

Nseg = size(TURB,1);
for i_s = 1:Nseg
    
    TURB(i_s).qv = TURB(i_s).qt - TURB(i_s).ql;
    
    [~,thv,th] = Tvthvth(TURB(i_s).thl,TURB(i_s).qv,TURB(i_s).ql,TURB(i_s).pbase);
    
    TURB(i_s).B  = detrend( g/mean(thv) *(thv-mean(thv)) );                % buoyancy
    TURB(i_s).Bt = detrend( g/mean(thv) *(th -mean(th)) );                 % pseudobuoyancy theta-part
    TURB(i_s).Bq = TURB(i_s).B - TURB(i_s).Bt;                             % pseudobuoyancy theta&q-part 
    
    TURB(i_s).cloud_mask = (TURB(i_s).ql>mask_ql_threshold);
    TURB(i_s).cloud_mask_ext = extend_mask(TURB(i_s).cloud_mask,mask_ext_factor,'proportional');
    TURB(i_s).cloud_fraction = sum(TURB(i_s).cloud_mask)/length(TURB(i_s).cloud_mask);
    TURB(i_s).cloud_fraction_ext = sum(TURB(i_s).cloud_mask_ext)/length(TURB(i_s).cloud_mask_ext);
    
    TURB(i_s).valid_fraction = 1;
    
end 


TURB = [TURB; TURB([TURB(:).level]=="cloud-base")];

for i_s = Nseg+1:length(TURB)
    TURB(i_s).level = "cloud-base-noclouds";
    TURB(i_s).valid_fraction = 1 - TURB(i_s).cloud_fraction_ext;
end



%% Compute structure functions for individual segments

disp('Structure functions ...')

% List of structure functions
sfc_vars = {'uuu',{'ux','ux','ux'};
            'vvu',{'vy','vy','ux'};
            'wwu',{'w','w','ux'};
            'wB', {'w','B'};
            'wBt',{'w','Bt'};
            'wBq',{'w','Bq'};
            'uu', {'ux','ux'};
            'vv', {'vy','vy'};
            'ww', {'w','w'};
            'wXuu',{'ux','ux','X_w'};
            'wXvv',{'vy','vy','X_w'};
            'wXww',{'w','w','X_w'};
            'wp', {'w','pi'}};


% Prepare data structures

S = rmfield(TURB,setdiff(fieldnames(TURB),{'level','dir','alt','valid_fraction','dr'})); % structure functions
V = S; % Variances of the form <dx^2>
L = S; % Integral length scales of increments
N = S; % Number of independent samples
U = S; % Uncertainties of structure functions


% Compute

Nseg = size(S,1);
Nvar = size(sfc_vars,1);
us = etd(clock,0,Nseg*Nvar,30);

for i_v = 1:Nvar
    
    var = sfc_vars{i_v,1};
    varlist = sfc_vars{i_v,2};
    
    midFlag = startsWith(varlist,'X_');  % factors of type w_X
    absFlag = startsWith(varlist,'A_');  % factors of type |dw|
    difFlag = ~(midFlag | absFlag);      % factors of type dw
    varlist(midFlag) = extractAfter(varlist(midFlag),'X_');
    varlist(absFlag) = extractAfter(varlist(absFlag),'A_');
    
    fprintf('S_%s = < ',var)
    if sum(midFlag), fprintf('%s_X ',  varlist{midFlag}), end
    if sum(absFlag), fprintf('|d %s| ',varlist{absFlag}), end
    if sum(difFlag), fprintf('d %s ',  varlist{difFlag}), end
    fprintf('>\n')

    for i_s = 1:Nseg
        
        % Data matrix
        A = cell2mat(cellfun(@(v) TURB(i_s).(v),varlist,'UniformOutput',false));
        
        % Mask cloudy points
        if S(i_s).level == "cloud-base-noclouds"
            A(logical(TURB(i_s).cloud_mask),:) = nan;
        end
        
        
        if ~any(all(isnan(A),1))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);
            S(i_s).N     = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                midI = 1; absI = 1; difI = 1;
                  
                % factors of type w_X
                if sum(midFlag)>0
                    midI = prod( 0.5*( circshift(A(:,midFlag),-r_lag) + A(:,midFlag) ) ,2);
                end
                
                % factors of type |dw|
                if sum(absFlag)>0
                    absI = prod( abs( circshift(A(:,absFlag),-r_lag) - A(:,absFlag) ) ,2);
                end
                
                % factors of type dw
                if sum(difFlag)>0
                    difI = prod( circshift(A(:,difFlag),-r_lag) - A(:,difFlag) ,2);
                end
                
                % increments
                I = midI.*absI.*difI; 
                
                % structure function
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        
                
                % number of non-nan samples
                S(i_s).N(r_lag) = sum(~isnan(I));            
                               
                % parameters needed for uncertainty
                if if_uncertainties
                    V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');      % variance of increments
                    L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr; % integral length scale for increments 
                end
                
            end
            
            % number of independent samples
%             N(i_s).(var) = S(i_s).N*S(i_s).dr./L(i_s).(var); 
            N(i_s).(var) = S(i_s).valid_fraction*seg_length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);

    end
    
end


%% Average at characteristic levels

disp('Average at characteristic levels ...')


avS = struct([]); avN = struct([]); avU = struct([]);

levels = unique([TURB(:).level],'stable');

Nlvl = numel(levels);
Nvar = size(sfc_vars,1);

for i_l = 1:Nlvl
    ind_l = find([S(:).level]'==levels(i_l));
    
    avS(i_l,1).level = levels(i_l);
    avS(i_l).number  = numel(ind_l);
    avS(i_l).alt     = mean([S(ind_l).alt]);
    avS(i_l).dr      = mean([S(ind_l).dr]);
    
    for i_v = 1:Nvar
        var = sfc_vars{i_v,1};
        
        avS(i_l,1).(var) = mean( vertcat(S(ind_l).(var)), 1, 'omitnan');
        avV              = mean( vertcat(V(ind_l).(var)), 1, 'omitnan');
        avN(i_l,1).(var) = sum(  vertcat(N(ind_l).(var)), 1);
        avU(i_l,1).(var) = sqrt( (avV - avS(i_l).(var).^2) ./ avN(i_l).(var) );
    end
    
    % s3
    if all(ismember({'uuu','vvu','wwu'},fieldnames(avS)))
        avS(i_l).S3 = 3*( avS(i_l).uuu + avS(i_l).vvu + avS(i_l).wwu );
        avU(i_l).S3 = 3*sqrt( avU(i_l).uuu.^2 + avU(i_l).vvu.^2 + avU(i_l).wwu.^2 );
    end
    
    % Transport integrand
    if all(ismember({'wXuu','wXvv','wXww'},fieldnames(avS)))
        avS(i_l).wXS2  = avS(i_l).wXuu + avS(i_l).wXvv + avS(i_l).wXww;
        avU(i_l).wXS2  = sqrt(avU(i_l).wXuu.^2 + avU(i_l).wXvv.^2 + avU(i_l).wXww.^2);
    end
    
end



%% Integrate / compensate

disp('Integrate or compensate ...')

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    avS(i_l).W  = 6*cumtrapz(r,avS(i_l).wB .*r.^2) ./ r.^3;
    avS(i_l).Wt = 6*cumtrapz(r,avS(i_l).wBt.*r.^2) ./ r.^3;
    avS(i_l).Wq = 6*cumtrapz(r,avS(i_l).wBq.*r.^2) ./ r.^3;
    avU(i_l).W  = 6*cumtrapz(r,avU(i_l).wB .*r.^2) ./ r.^3;
    avU(i_l).Wt = 6*cumtrapz(r,avU(i_l).wBt.*r.^2) ./ r.^3;
    avU(i_l).Wq = 6*cumtrapz(r,avU(i_l).wBq.*r.^2) ./ r.^3;
    
    avS(i_l).Tiu = 3*cumtrapz(r,avS(i_l).wXS2.*r.^2) ./ r.^3;
    avU(i_l).Tiu = 3*cumtrapz(r,avU(i_l).wXS2.*r.^2) ./ r.^3;
    
    avS(i_l).Tip = 6*cumtrapz(r,avS(i_l).wp.*r.^2) ./ r.^3;
    avU(i_l).Tip = 6*cumtrapz(r,avU(i_l).wp.*r.^2) ./ r.^3;
    
    avS(i_l).uuu3r = 3*avS(i_l).uuu./r;
    avS(i_l).vvu3r = 3*avS(i_l).vvu./r;
    avS(i_l).wwu3r = 3*avS(i_l).wwu./r;
    avU(i_l).uuu3r = 3*avU(i_l).uuu./r;
    avU(i_l).vvu3r = 3*avU(i_l).vvu./r;
    avU(i_l).wwu3r = 3*avU(i_l).wwu./r;
    
    avS(i_l).S3r = avS(i_l).S3./r;
    avU(i_l).S3r = avU(i_l).S3./r;
    
    avS(i_l).mS3 = -avS(i_l).S3;
    avU(i_l).mS3 = avU(i_l).S3;
    
    avS(i_l).WmS3r = avS(i_l).W - avS(i_l).S3r;
    avU(i_l).WmS3r = avU(i_l).W + avU(i_l).S3r;
    
    avS(i_l).WrmS3 = avS(i_l).W.*r - avS(i_l).S3;
    avU(i_l).WrmS3 = avU(i_l).W.*r + avU(i_l).S3;
end




%% Estimate transport (analytical method)

disp('Estimate transport ...')

interp_method = 'pchip';

level_sets = {'cloud-base','top-subcloud','mid-subcloud','near-surface';
     'cloud-base-noclouds','top-subcloud','mid-subcloud','near-surface'};

for i_ls = 1:size(level_sets,1) 
    [~,ind_l] = ismember(level_sets(i_ls,:),[avS(:).level]);

    alt_vec   = [avS(ind_l).alt];

    [alt_vec,ind_sort] = sort(alt_vec);
    ind_l = ind_l(ind_sort);

    
    Tiu_interp = interp1(alt_vec,vertcat(avS(ind_l).Tiu),interp_method,'pp');       
    if if_uncertainties
        u_Tiu_interp = interp1(alt_vec,vertcat(avU(ind_l).Tiu),interp_method,'pp');
    else
        u_Tiu_interp.coefs = nan(size(Tiu_interp.coefs));
        u_Tiu_interp.dim = Tiu_interp.dim;
    end

    Nll = numel(ind_l);
    for ii_l = 1:Nll-1
        i_l = ind_l(ii_l);
        avS(i_l).Tu =   Tiu_interp.coefs( 1+Tiu_interp.dim*(ii_l-1) : Tiu_interp.dim*ii_l, 3 )';
        avU(i_l).Tu = u_Tiu_interp.coefs( 1+u_Tiu_interp.dim*(ii_l-1) : u_Tiu_interp.dim*ii_l, 3 )';
    end
    i_l = ind_l(end);        
    dz = Tiu_interp.breaks(end) - Tiu_interp.breaks(end-1);
    CFS = Tiu_interp.coefs( 1 + Tiu_interp.dim*(Nll-2) : end, : );
    avS(i_l).Tu = ( 3*CFS(:,1)*dz.^2 + 2*CFS(:,2)*dz + CFS(:,3) )';
    uCFS = u_Tiu_interp.coefs( 1 + u_Tiu_interp.dim*(Nll-2) : end, : );
    avU(i_l).Tu = ( 3*uCFS(:,1)*dz.^2 + 2*uCFS(:,2)*dz + uCFS(:,3) )';
    
    
    Tip_interp = interp1(alt_vec,vertcat(avS(ind_l).Tip),interp_method,'pp');       
    if if_uncertainties
        u_Tip_interp = interp1(alt_vec,vertcat(avU(ind_l).Tip),interp_method,'pp');
    else
        u_Tip_interp.coefs = nan(size(Tip_interp.coefs));
        u_Tip_interp.dim = Tip_interp.dim;
    end

    Nll = numel(ind_l);
    for ii_l = 1:Nll-1
        i_l = ind_l(ii_l);
        avS(i_l).Tp =   Tip_interp.coefs( 1+Tip_interp.dim*(ii_l-1) : Tip_interp.dim*ii_l, 3 )';
        avU(i_l).Tp = u_Tip_interp.coefs( 1+u_Tip_interp.dim*(ii_l-1) : u_Tip_interp.dim*ii_l, 3 )';
    end
    i_l = ind_l(end);        
    dz = Tip_interp.breaks(end) - Tip_interp.breaks(end-1);

    CFS = Tip_interp.coefs( 1 + Tip_interp.dim*(Nll-2) : end, : );
    avS(i_l).Tp = ( 3*CFS(:,1)*dz.^2 + 2*CFS(:,2)*dz + CFS(:,3) )';
    uCFS = u_Tip_interp.coefs( 1 + u_Tip_interp.dim*(Nll-2) : end, : );
    avU(i_l).Tp = ( 3*uCFS(:,1)*dz.^2 + 2*uCFS(:,2)*dz + uCFS(:,3) )';
end


% Dependent parameters

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    avS(i_l).T = avS(i_l).Tu + avS(i_l).Tp;
    avU(i_l).T = sqrt( avU(i_l).Tu.^2 + avU(i_l).Tp.^2 );
    
    avS(i_l).WmS3rmT = avS(i_l).W - avS(i_l).S3r - avS(i_l).T;
    avU(i_l).WmS3rmT = avU(i_l).W + avU(i_l).S3r + avU(i_l).T;
    
    avS(i_l).WrmS3mTr = avS(i_l).W.*r - avS(i_l).S3 - avS(i_l).T.*r;
    avU(i_l).WrmS3mTr = avU(i_l).W.*r + avU(i_l).S3 + avU(i_l).T.*r;
end



%% Estimate dissipation rate

disp('Estimate dissipation rate ...')

edr_vars = {'uu','vv','ww','mS3','WrmS3','WrmS3mTr'};
slps = [2/3 2/3 2/3 1 1 1];
cons = [2.0 2.6 2.6 4.0 4.0 4.0];

Nvar = numel(edr_vars);

for i_v = 1:Nvar
    var = edr_vars{i_v};
    [avS,avU] = calc_edr(avS,avU,edr_vars{i_v},edr_fit_range,cons(i_v),...
        'Slope',slps(i_v),'Method','direct','dC',cons(i_v)*0.15);
end


% Compensate for plotting

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    for i_v = 1:Nvar
        var = edr_vars{i_v};
        avS(i_l).([var,'_c']) = avS(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
        avU(i_l).([var,'_c']) = avU(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
    end
    
    avS(i_l).edr_S2_4 = avS(i_l).edr_S2*4;
    avU(i_l).edr_S2_4 = avU(i_l).edr_S2*4;
end




%% Save/load

save('S_bomex.mat','r_max','edr_fit_range','levels','mask*','sfc_vars',...
     'S','U','L','N','V','avS','avU','avN',...
    'plotpath')


% load('S_eureca_1km.mat')




%% PLOTS

xlim = [dr r_max];
ylim = [1e-6 5e-3];
% ylim = [-inf inf];

Npoints = 40;

Nlvl = size(avS,1);


% S2

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim,'YLim',[5e-4 5e-3]);
    if i_l>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}C^{-1} $','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s2_',levels{i_l}],'-dpng','-r300')
end
    

%% W

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','Wt','Wq'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W$','$W_\theta$','$W_q$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'W_',levels{i_l}],'-dpng','-r300')
end


%% S3

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'S3r'},[1 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$S_3 r^{-1}$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3r_',levels{i_l}],'-dpng','-r300')
end


%% T

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'T','Tu','Tp'},[4 5 9],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$T$','$T_u$','$T_p$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'T_',levels{i_l}],'-dpng','-r300')
end


%% Balance: W, S3/r, T, epsilon

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','-Tu','-T','edr_S2_4'},[7 1 5 4 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$-T$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_l}],'-dpng','-r300')
end


%% Total sum: W - S3/r - T

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3rmT','edr_S2_4'},[2 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W-S_3 r^{-1}-T$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
    print(fig,[plotpath,filesep,'total_',levels{i_l}],'-dpng','-r300')
end


%% S3 components

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uuu3r','vvu3r','wwu3r'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'uuu','vvu','wwu'},'Interpreter','latex','Location','best')
    end
    ylabel('$3S_3 r^{-1}\,[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3_comp_',levels{i_l}],'-dpng','-r300')
end


%% Fitting dissipation

for i_l = 1:Nlvl
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim,'YLim',[1e-3 2e-2]);
    if i_l>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2fit_',levels{i_l}],'-dpng','-r300')
end


