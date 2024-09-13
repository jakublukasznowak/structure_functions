
%% Settings

% Cloud mask

mask_rh_threshold = 98;   % Relative humidity threshold for cloud mask
mask_ext_factor   = 1;    % Mask extension factor to account for cloud shells
mask_nan_replace  = true; % Assume cloud if masks are nan


% Correlation functions

r_max = 10e3; % [m]



%% Load list

% List of levels
levels  = {'cloud-base','top-subcloud','mid-subcloud','near-surface','cloud-base-noclouds'};

% List of flights
flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of variables from turbulent moments dataset
mom_vars = {'alt','time_start','time_end',...
    'MEAN_P','MEAN_THETA','MEAN_MR','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset
turb_vars = {'time','UX_DET','VY_DET','W_DET','T_DET','MR_DET'};

% List of variables from cloud microphysics dataset
mcph_vars = {'time','LWC','CLOUD_mask'};


% List of segments to exclude
% exclude_seg = ["RF19","S1";   % (surface) altitude changes + temperature drop
%                "RF09","L1B";  % (mid-subcloud) temperature drop + wind direction variations + TAS changes
%                "RF11","L2B";  % (mid-subcloud) rain + TAS changes
%                "RF17","L2A";  % (mid-subcloud) heavy rain + temperature drop + TAS changes
%                "RF11","L1A";  % (top-subcloud) heavy rain
%                "RF17","L1B";  % (top-subcloud) heavy rain + temperature variations + TAS changes
%                "RF17","R2B";  % (cloud-base) heavy rain
%                "RF18","R1B"]; % (cloud-base) heavy rain
exclude_seg = [];


% List of segment orientations
dirs2 = {'cross','along'};



%% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures',filesep,'correlations'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Load datasets

% Flight segmentation

disp('Load flight segmentation ...')
SEG = load_atr_seg(datapath,'v1.9','longlegs');                          
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);


% Exclude problematic segments

Nexc = size(exclude_seg,1);
for i_e = 1:Nexc
    ind_s = find( SEG.flight==exclude_seg(i_e,1) & SEG.name==exclude_seg(i_e,2) );
    SEG(ind_s,:) = [];
end


% Mean values and moments

disp('Load mean values and moments ...')
[MOM,mom_info] = load_atr_mom(datapath,'L3','v1.9','longlegs',mom_vars); 
MOM = join(SEG,MOM,'Keys',{'start','end'});


% Turbulent fluctuations

disp('Load turbulence data ...')
[TURB,turb_info] = load_atr_turb(MOM,datapath,'L3','v1.9',turb_vars);

MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


% Microphysics

disp('Load microphysics data ...')
MCPH = load_atr_pma(MOM,datapath,mcph_vars);


Nseg = size(MOM,1);
fprintf('Number of segments loaded: %d\n',Nseg)



%% Compute auxilliary parameters

disp('Compute thermodynamics and cloud mask ...')


% Thermodynamics (including buoyancy)

TURB = calc_thermo(TURB,MOM);


% Cloud masks

[TURB,MOM] = cloud_mask(TURB,MOM,MCPH,mask_rh_threshold,mask_ext_factor,mask_nan_replace);

for i_s = 1:Nseg
    TURB(i_s).OR_mask(isnan(TURB(i_s).OR_mask)) = mask_nan_replace;
end
MOM.valid_fraction = 1 - MOM.OR_cloud_fraction;


% Duplicate cloud-base segments

MOM = [MOM; MOM(MOM.level=="cloud-base",:)];
MOM.level(Nseg+1:end) = "cloud-base-noclouds";
MOM.valid_fraction(MOM.level~="cloud-base-noclouds") = 1;

TURB = [TURB; TURB([TURB(:).level]=="cloud-base")];


% Assume fixed distance between points

dr = 4;
MOM.dr(:) = dr;
r_maxlag = r_max/dr;



%% Compute correlation functions for individual segments

disp('Compute correlation functions ...')

% List of structure functions
corr_vars = {'uu',{'UX_DET','UX_DET'};
             'vv',{'VY_DET','VY_DET'};
             'ww',{'W_DET','W_DET'}};
         

% Prepare data structures

R = table2struct(MOM(:,{'flight','name','level','dir2','alt','dr','length','valid_fraction'})); % correlation functions
L = R; % Integral length scales


% Compute

Nvar = size(corr_vars,1);
Nseg = size(R,1);
us = etd(clock,0,Nseg*Nvar,30);

for i_v = 1:Nvar
    var = corr_vars{i_v,1};

    for i_s = 1:Nseg
        dr = MOM.dr(i_s);
        
        % Data vectors
        x = TURB(i_s).(corr_vars{i_v,2}{1});
        y = TURB(i_s).(corr_vars{i_v,2}{2});
        
        % Mask cloudy points
        if R(i_s).level == "cloud-base-noclouds"
            x(logical(TURB(i_s).OR_mask)) = nan;
            y(logical(TURB(i_s).OR_mask)) = nan;
        end
        
        % Correlation function
        if ~any(all(isnan([x y]),1))
            [LS,rho] = int_ls_short(x,y,'Method','e-decay','MaxLag',r_maxlag,'dr',dr);
            R(i_s).(var) = rho';
            L(i_s).(var) = LS;
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);
    end
    
end



%% Average at characteristic levels

disp('Average at characteristic levels ...')

avR = struct([]); avL = struct([]); lvL = struct([]);

Nlvl = numel(levels);
Ndir = numel(dirs2);
Nvar = size(corr_vars,1);

for i_l = 1:Nlvl
    for i_d = 1:Ndir
        ind_ld = find([R(:).level]'==levels{i_l} & [R(:).dir2]'==dirs2{i_d});
        
        avR(i_l,i_d).level = string(levels{i_l});
        avR(i_l,i_d).dir2  = string(dirs2{i_d});
        avR(i_l,i_d).number  = numel(ind_ld);
        avR(i_l,i_d).alt     = mean([R(ind_ld).alt]);
        avR(i_l,i_d).dr      = mean([R(ind_ld).dr]);
        avR(i_l,i_d).length  = mean([R(ind_ld).length]);
    
        if ~isempty(ind_ld)
            for i_v = 1:Nvar
                var = corr_vars{i_v,1};

                avR(i_l,i_d).(var) = mean( vertcat(R(ind_ld).(var)), 1);
                avL(i_l,i_d).(var) = mean( vertcat(L(ind_ld).(var)), 1);
                lvL(i_l,i_d).(var) = int_ls_long(avR(i_l,i_d).(var),'Method','e-decay','dr',avR(i_l,i_d).dr);
            end
        end
    end
end



%% Save/load

save('R_eureca.mat','r_max','levels','dirs2','exclude_seg','mask*','corr_vars',...
    'MOM','R','L','avR','avL','lvL',...
    'plotpath')

% load('R_eureca.mat')



%% PLOTS

lw = 1;
mks = 10;

[f,~,co,ls,mk] = fig16x12;
close(f)
co = repmat(co,100,1);

Nlvl = numel(levels);
Ndir = numel(dirs2);
Nvar = size(corr_vars,1);
Nseg = size(R,1);


%% R segments + levels

for i_v = 1:Nvar
    var = corr_vars{i_v,1};
    
    for i_l = 1:Nlvl
        lvl = levels{i_l};
    
        fig = figure('Units','normalized','Position',[0 0.5 1 0.4]);
        grid on, hold on
        i_c = 1; legstr = [];
        
        for i_s = 1:Nseg
            dr = R(i_s).dr;
            rv = (0:length(R(i_s).(var))-1)*dr/1000;
            
            if R(i_s).level==lvl
                c = co(i_c,:);
                ind_d = find(R(i_s).dir2==dirs2);
                plot(rv,R(i_s).(var),ls{ind_d},'Color',c,'LineWidth',lw);
                plot(L(i_s).(var)/1000,...
                    interp1(rv,R(i_s).(var),L(i_s).(var)/1000,'linear','extrap'),...
                    mk{ind_d},'Color',c,'MarkerSize',mks-2,'HandleVisibility','off');
                i_c = i_c + 1;
                legstr = vertcat(legstr,R(i_s).flight+" "+R(i_s).name);
            end
        end
        
        i_c = 8; c = co(i_c,:);
        dr = avR(i_l).dr;
        rv = (0:length(R(i_s).(var))-1)*dr/1000;
        for i_d = 1:Ndir
            avr = avR(i_l,i_d).(var);
            if ~isempty(avr)
                plot(rv,avr,ls{i_d},'Color',c,'LineWidth',lw+2)
                plot(avL(i_l,i_d).(var)/1000,...
                    interp1(rv,avr,avL(i_l,i_d).(var)/1000,'linear','extrap'),...
                    mk{i_d},'Color',c,...
                    'MarkerSize',mks+2,'HandleVisibility','off');
                plot(lvL(i_l,i_d).(var)/1000,...
                    interp1(rv,avr,lvL(i_l,i_d).(var)/1000,'linear','extrap'),...
                    mk{i_d},'Color',c,'MarkerFaceColor',c,...
                    'MarkerSize',mks+2,'HandleVisibility','off');
                legstr = vertcat(legstr,"av "+avR(i_l,i_d).dir2);
            end
        end
        
        plot([0,rv(end)],exp(-1)*[1 1],'LineStyle','--','LineWidth',lw,'Color','black')
        plot([0,rv(end)],[0 0],'LineStyle','--','LineWidth',lw,'Color','black')
        
        set(gca,'XLim',[0 r_max/1000])
        xlabel('r [km]')
        ylabel(['R ',var])
        legend(legstr,'Location','eastoutside')
        title([lvl,' (solid: cross, dashed: along)'])
        print(fig,[plotpath,filesep,'R_',lower(var),'_',levels{i_l}],'-dpng','-r300')
    end
    
end
 

%% R levels

for i_v = 1:Nvar
    var = corr_vars{i_v,1};

    fig = figure('Units','normalized','Position',[0 0 1 0.4]);
    grid on, hold on
    
    for i_d = 1:Ndir
        for i_l = 1:Nlvl
            
            if ~isempty(avR(i_l,i_d).(var))
                Rv = avR(i_l,i_d).(var);
                rv = (0:length(Rv)-1)*avR(i_l,i_d).dr/1000;
                plot(rv,Rv,ls{i_d},'Color',co(i_l,:),'LineWidth',lw,'HandleVisibility','off')
            end            
        end
    end
    for i_l = 1:Nlvl, plot(nan,nan,ls{1},'Color',co(i_l,:)), end
    for i_d = 1:Ndir, plot(nan,nan,ls{i_d},'Color','black'), end
    
    plot([0,rv(end)],exp(-1)*[1 1],'--','LineWidth',lw,'Color','black')
    plot([0,rv(end)],[0 0],        '--','LineWidth',lw,'Color','black')
    
    set(gca,'XLim',[0 r_max/1000])
    xlabel('r [km]')
    ylabel(['R ',var,' [m]'])
    legend(horzcat(levels,dirs2),'Location','eastoutside')
    print(fig,[plotpath,filesep,'avR_',lower(var)],'-dpng','-r300')
end


%% L levels

for i_v = 1:Nvar
    var = corr_vars{i_v,1};
    
    fig = figure('Units','normalized','Position',[0 0 0.4 0.5]);
    grid on, hold on
    
    for i_d = 1:Ndir
        for i_l = 1:Nlvl
            if ~isempty(avL(i_l,i_d).(var))
                plot(avL(i_l,i_d).(var),avR(i_l,i_d).alt,mk{i_d},...
                    'Color',co(i_l,:),'MarkerFaceColor',co(i_l,:),...
                    'MarkerSize',mks,'HandleVisibility','off')
            end
            if ~isempty(lvL(i_l,i_d).(var))
                plot(lvL(i_l,i_d).(var),avR(i_l,i_d).alt,mk{i_d},...
                    'Color',co(i_l,:),'MarkerSize',mks,'HandleVisibility','off')
            end
        end
    end
        
    for i_l=1:Nlvl, plot(nan,nan,mk{1},'Color',co(i_l,:),'MarkerSize',mks), end
    for i_d=1:Ndir, plot(nan,nan,mk{i_d},'Color','black','MarkerSize',mks), end
    plot(nan,nan,mk{1},'Color','black','MarkerFaceColor','black','MarkerSize',mks)
    plot(nan,nan,mk{1},'Color','black','MarkerSize',mks)
    
    if ismember(var,{'uu','vv'})
        set(gca,'XLim',[0 2500])
    end
    ylabel('z [m]')
    xlabel(['L ',var,' [m]'])
    legend(horzcat(levels,dirs2,{'av seg-L','level-L'}),'Location','northwest')
    print(fig,[plotpath,filesep,'avL_',lower(var)],'-dpng','-r300')
end

