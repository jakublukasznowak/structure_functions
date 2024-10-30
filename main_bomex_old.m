
tq = 6*3600; % [s]
zq = [60 250 450 600]; % m
levels = {'near-surface','mid-subcloud','top-subcloud','cloud-base'};
ind_x = 10:10:380;

r_max = 1000;
edr_fit_range = [30 300];
if_uncertainties = true;
transport_method = 'interp';

outputfile = 'S_bomex_1km';



%% Load

vars = {'qt','thl','u','v','w'};

Nvar = numel(vars);
Nlvl = numel(zq);
Nxq = numel(ind_x);


% Paths

datapath = [myprojectpath,filesep,'data_bomex_2406',filesep];
addpath(genpath(myprojectpath))

plotpath = [myprojectpath,filesep,'figures_LES'];
if ~isfolder(plotpath), mkdir(plotpath), end


% Coordinates

t = ncread([datapath,'qt.nc'],'time');
x = ncread([datapath,'qt.nc'],'x');
y = ncread([datapath,'qt.nc'],'y');
z = ncread([datapath,'qt.nc'],'z');
xh = ncread([datapath,'u.nc'],'xh');
yh = ncread([datapath,'v.nc'],'yh');
zh = ncread([datapath,'w.nc'],'zh');

% Slices in time
ind_t = find(t>=tq,1,'first');

% Slices in altitude
ind_z = zeros(1,Nlvl);
for i_z = 1:Nlvl
    ind_z(i_z) = find(z>=zq(i_z),1,'first');
end
zqr = z(ind_z);

% Slices in x
xq = x(ind_x);


% Profiles

disp('Load profiles ...')

for i_v = 1:Nvar
    var = vars{i_v};
    prof.(var) = squeeze(mean(ncread([datapath,var,'.nc'],var),[1 2]));
end

% Plot profiles at tq
figure, hold on, grid on
set(gca,'XColor','blue')
plot(prof.qt(:,ind_t)*1e3,z,'b')
xlabel('qt [g/kg]'), ylabel('z [m]')
axes('Parent',gcf,'Position',get(gca,'Position'),'Color','none',...
    'XAxisLocation','top','YAxisLocation','right','XColor','red'); hold on
plot(prof.thl(:,ind_t),z,'r')
xlabel('thl [deg K]'), ylabel('z [m]')
for i_z = 1:Nlvl
    plot(get(gca,'XLim'),zqr(i_z)*[1 1],'g','LineWidth',1)
end
print(gcf,[plotpath,filesep,'profile'],'-r300','-dpng')


% Horizontal planes

disp('Load horizontal planes ...')

plane = struct([]);

for i_z = 1:Nlvl
    for i_v = 1:Nvar
        var = vars{i_v};
        if strcmp(var,'w'), zn = 2;
        else, zn = 1; end
        plane(i_z,1).(var) = ncread([datapath,var,'.nc'],var,[1 1 ind_z(i_z) ind_t],[Inf Inf zn 1]);
    end
end

% Interpolate to basic grid
for i_z = 1:Nlvl
    plane(i_z).u = interp1(xh,plane(i_z).u,x);
    plane(i_z).v = shiftdim(interp1(yh,shiftdim(plane(i_z).v,1),y),1);
    plane(i_z).w = shiftdim(interp1(zh([ind_z(i_z) ind_z(i_z)+1]),shiftdim(plane(i_z).w,2),z(ind_z(i_z))),1);
end

% Remove nan-lines at edges
for i_z = 1:Nlvl
    for i_v = 1:Nvar
        var = vars{i_v};
        plane(i_z).(var) = plane(i_z).(var)(1:end-1,1:end-1);
        x = x(1:end-1); y = y(1:end-1);
    end
end

% Plot planes
var = 'u';
for i_z = 1:Nlvl
    figure, hold on
    imagesc(x,y,plane(i_z).(var)')
    xlim([0 x(end)]), ylim([0 y(end)])
    colormap hot
    colorbar
    caxis([-0.5 3.5])
    for i_x = 1:Nxq
        plot(xq(i_x)*[1 1],get(gca,'YLim'),'g','LineWidth',1)
    end
    title(sprintf('%s at z=%.0fm t=%.0fs',var,zqr(i_z),tq))
    print(gcf,[plotpath,filesep,'plane_',var,'_',levels{i_z}],'-r300','-dpng')
end


% Cut out segments

disp('Cut out segments ...')

TURB = cell(Nlvl,1);

for i_z = 1:Nlvl
    TURB{i_z} = struct([]);
 
    for i_x = 1:Nxq
        TURB{i_z}(i_x,1).level = string(levels{i_z});
        TURB{i_z}(i_x).ind_z = ind_z(i_z);
        TURB{i_z}(i_x).z = zqr(i_z);
        TURB{i_z}(i_x).alt = zqr(i_z);
        TURB{i_z}(i_x).ind_x = ind_x(i_x);
        TURB{i_z}(i_x).x = xq(i_x);
        TURB{i_z}(i_x).length = y(end);
        
        for i_v = 1:Nvar
            var = vars{i_v};
            TURB{i_z}(i_x).(var) = plane(i_z).(var)(ind_x(i_x),:)';
        end 
    end
end

TURB = vertcat(TURB{:});



%% Thermodynamics

disp('Compute thermodynamics ...')

% Constants
Mv = 18.015;     % [g/mol]
Md = 28.9647;
R  = 8.31432;    % [J/mol/K]
Rd = R/Md*1000;  % [J/kg/K]
Rv = R/Mv*1000;
ep = Rd/Rv;
cpd= 1003.6;     % [J/kg/K]
Lv = 2.5008e6;   % [J/kg]
T0 = 273.15;     % [K]
g  = 9.81;       % [m/s2]

% Approximation of Clausius-Clapeyron
es = @(t) 6.112*exp(17.67*t./(t+243.5));  % hPa


Nseg = size(TURB,1);

for i_s = 1:Nseg
    qv = TURB(i_s).qt; % -ql
%     ql = TURB(i_s).ql;
    th = TURB(i_s).thl; % +Lv/cpd*ql
    thv = th.*(1+(1/ep-1)*qv); % -ql
    
    TURB(i_s).B  = detrend( g/mean(thv) *(thv-mean(thv)) );                % buoyancy
    TURB(i_s).Bt = detrend( g/mean(thv) *(th -mean(th)) );                 % pseudobuoyancy theta-part
    TURB(i_s).Bq = detrend( g/mean(thv) *(1/ep-1)*(th.*qv-mean(th.*qv)) ); % -th*ql-mean(th*ql) % pseudobuoyancy theta&q-part 
end


for i_s = 1:Nseg
    TURB(i_s).u = detrend(TURB(i_s).u);
    TURB(i_s).v = detrend(TURB(i_s).v);
    TURB(i_s).w = detrend(TURB(i_s).w);
end



%% Structure functions

disp('Compute structure functions ...')

dr = unique(diff(yh));
r_maxlag = r_max/dr;

% List of structure functions
sfc_vars = {'uuu',{'v','v','v'};
            'vvu',{'u','u','v'};
            'wwu',{'w','w','v'};
            'wB', {'w','B'};
            'wBt',{'w','Bt'};
            'wBq',{'w','Bq'};
            'uu', {'v','v'};
            'vv', {'u','u'};
            'ww', {'w','w'};
            'wXuu',{'v','v','X_w'};
            'wXvv',{'u','u','X_w'};
            'wXww',{'w','w','X_w'}};
 
% Prepare data structures
S = rmfield(TURB,setdiff(fieldnames(TURB),{'level','dr','x','z','alt','ind_x','ind_z','length'}));
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
        
        TURB(i_s).dr = dr;
        S(i_s).dr = TURB(i_s).dr;
        S(i_s).valid_fraction = 1;
        
        % Data matrix
        A = cell2mat(cellfun(@(v) TURB(i_s).(v),varlist,'UniformOutput',false));
        
        % Mask cloudy points
        if S(i_s).level == "cloud-base-noclouds"
            A(logical(TURB(i_s).OR_mask),:) = nan;
        end
        
        Lt = length(TURB(i_s).u);
        
        if ~any(all(isnan(A),1))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);
            S(i_s).N     = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                midI = 1; absI = 1; difI = 1;
                  
                % factors of type w_X
                if sum(midFlag)>0
                    midI = prod( 0.5*(A(r_lag+1:Lt,midFlag)+A(1:Lt-r_lag,midFlag)) ,2);
                end
                
                % factors of type |dw|
                if sum(absFlag)>0
                    absI = prod( abs( A(r_lag+1:Lt,absFlag)-A(1:Lt-r_lag,absFlag) ) ,2);
                end
                
                % factors of type dw
                if sum(difFlag)>0
                    difI = prod( A(r_lag+1:Lt,difFlag)-A(1:Lt-r_lag,difFlag) ,2);
                end
                
                % increments
                I = midI.*absI.*difI; 
                
                % structure function
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        
                
                % parameters needed for uncertainty
                if if_uncertainties
                    V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');      % variance of increments
                    L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr; % integral length scale for increments
                    S(i_s).N    (r_lag) = sum(~isnan(I));            % number of valid samples
                end
                
            end
            
            % number of independent samples
%             N(i_s).(var) = S(i_s).N*dr./L(i_s).(var); 
            N(i_s).(var) = S(i_s).valid_fraction*S(i_s).length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);

    end
    
end



%% Average at characteristic levels

disp('Average at characteristic levels ...')


avS = struct([]); avN = struct([]); avU = struct([]);

Nlvl = numel(levels);
Nvar = size(sfc_vars,1);

for i_z = 1:Nlvl
    ind_l = find([S(:).level]'==levels{i_z});
    
    avS(i_z,1).level = string(levels{i_z});
    avS(i_z).number  = numel(ind_l);
    avS(i_z).z       = mean([S(ind_l).z]);
    avS(i_z).x       = mean([S(ind_l).x]);
    avS(i_z).dr      = mean([S(ind_l).dr]);
    avS(i_z).alt     = mean([S(ind_l).alt]);
    
    for i_v = 1:Nvar
        var = sfc_vars{i_v,1};
        
        avS(i_z,1).(var) = mean( vertcat(S(ind_l).(var)), 1);
        avV              = mean( vertcat(V(ind_l).(var)), 1);
%         avN(i_z,1).(var) = sum( vertcat(N(ind_l).(var)) ./ vertcat(L(ind_l).(var)) );
        avN(i_z,1).(var) = sum( repmat([S(ind_l).valid_fraction]'.*[S(ind_l).length]',1,r_maxlag) ...
            ./ vertcat(L(ind_l).(var)) );
        avU(i_z,1).(var) = sqrt( (avV - avS(i_z).(var).^2) ./ avN(i_z).(var) );
    end
    
    % s3
    if all(ismember({'uuu','vvu','wwu'},fieldnames(avS)))
        avS(i_z).S3 = 3*( avS(i_z).uuu + avS(i_z).vvu + avS(i_z).wwu );
        avU(i_z).S3 = 3*sqrt( avU(i_z).uuu.^2 + avU(i_z).vvu.^2 + avU(i_z).wwu.^2 );
    end
    
    % Transport integrand
    if all(ismember({'wXuu','wXvv','wXww'},fieldnames(avS)))
        avS(i_z).wXS2  = avS(i_z).wXuu + avS(i_z).wXvv + avS(i_z).wXww;
        avU(i_z).wXS2  = sqrt(avU(i_z).wXuu.^2 + avU(i_z).wXvv.^2 + avU(i_z).wXww.^2);
    end
    
    % Consistency check for wB = wBt + wBq
    if all(ismember({'wBt','wBq'},fieldnames(avS)))
        avS(i_z).wQ  = avS(i_z).wBt + avS(i_z).wBq;
        avU(i_z).wQ  = sqrt(avU(i_z).wBt.^2 + avU(i_z).wBq.^2);
    end
    
end



%% Integrate / compensate

disp('Integrate or compensate ...')

Nlvl = size(avS,1);

for i_z = 1:Nlvl
    r = (1:length(avS(i_z).uuu))*avS(i_z).dr;
    
    avS(i_z).W  = 6*cumtrapz(r,avS(i_z).wB .*r.^2) ./ r.^3;
    avS(i_z).Wt = 6*cumtrapz(r,avS(i_z).wBt.*r.^2) ./ r.^3;
    avS(i_z).Wq = 6*cumtrapz(r,avS(i_z).wBq.*r.^2) ./ r.^3;
    avU(i_z).W  = 6*cumtrapz(r,avU(i_z).wB .*r.^2) ./ r.^3;
    avU(i_z).Wt = 6*cumtrapz(r,avU(i_z).wBt.*r.^2) ./ r.^3;
    avU(i_z).Wq = 6*cumtrapz(r,avU(i_z).wBq.*r.^2) ./ r.^3;
    
    avS(i_z).Ti = 3*cumtrapz(r,avS(i_z).wXS2.*r.^2) ./ r.^3;
    avU(i_z).Ti = 3*cumtrapz(r,avU(i_z).wXS2.*r.^2) ./ r.^3;
    
    avS(i_z).uuu3r = 3*avS(i_z).uuu./r;
    avS(i_z).vvu3r = 3*avS(i_z).vvu./r;
    avS(i_z).wwu3r = 3*avS(i_z).wwu./r;
    avU(i_z).uuu3r = 3*avU(i_z).uuu./r;
    avU(i_z).vvu3r = 3*avU(i_z).vvu./r;
    avU(i_z).wwu3r = 3*avU(i_z).wwu./r;
    
    avS(i_z).S3r = avS(i_z).S3./r;
    avU(i_z).S3r = avU(i_z).S3./r;
    
    avS(i_z).mS3 = -avS(i_z).S3;
    avU(i_z).mS3 = avU(i_z).S3;
    
    avS(i_z).WmS3r = avS(i_z).W - avS(i_z).S3r;
    avU(i_z).WmS3r = avU(i_z).W + avU(i_z).S3r;
    
    avS(i_z).WrmS3 = avS(i_z).W.*r - avS(i_z).S3;
    avU(i_z).WrmS3 = avU(i_z).W.*r + avU(i_z).S3;
end



%% Estimate transport

disp('Estimate transport ...')


if strcmp(transport_method,'diff') % from differences between levels
    
    level_pairs = {'cloud-base','top-subcloud';
                   'top-subcloud','mid-subcloud';
                   'mid-subcloud','near-surface'};

    for i_lp = 1:size(level_pairs,1)
        i_l1 = find([avS(:).level]==level_pairs{i_lp,1});
        i_l2 = find([avS(:).level]==level_pairs{i_lp,2});

        avS(i_l1).T = (avS(i_l1).Ti-avS(i_l2).Ti) / (avS(i_l1).alt-avS(i_l2).alt);
        avU(i_l1).T = sqrt(avU(i_l1).Ti.^2+avU(i_l2).Ti.^2) / (avS(i_l1).alt-avS(i_l2).alt);
    end

    
elseif strcmp(transport_method,'interp') % from interpolated smooth function

    alt_increment = 25; % m
    interp_method = 'pchip';

    level_sets = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};
    
    for i_ls = 1:size(level_sets,1) 
        [~,ind_l] = ismember(level_sets(i_ls,:),[avS(:).level]);
        
        alt_vec   = [avS(ind_l).z];
        alt_query = reshape( [alt_vec-alt_increment; alt_vec+alt_increment] ,1,[]);
        
        Ti_interp = interp1(alt_vec,vertcat(avS(ind_l).Ti),alt_query,interp_method);
        if if_uncertainties
            u_Ti_interp = interp1(alt_vec,cat(1,avU(ind_l).Ti),alt_query,interp_method);
        else
            u_Ti_interp = nan(size(Ti_interp));
        end

        for ii_l = 1:numel(ind_l)
            i_z = ind_l(ii_l);
            avS(i_z).T = ( Ti_interp(2*ii_l,:) - Ti_interp(2*ii_l-1,:) )/alt_increment/2;
            avU(i_z).T = sqrt( u_Ti_interp(2*ii_l,:).^2 + u_Ti_interp(2*ii_l-1,:).^2 )/alt_increment/2;
        end
    end
    
end


% Dependent parameters

Nlvl = size(avS,1);

for i_z = 1:Nlvl
    r = (1:length(avS(i_z).uuu))*avS(i_z).dr;

    avS(i_z).WmS3rmT = avS(i_z).W - avS(i_z).S3r - avS(i_z).T;
    avU(i_z).WmS3rmT = avU(i_z).W + avU(i_z).S3r + avU(i_z).T;
    
    avS(i_z).WrmS3mTr = avS(i_z).W.*r - avS(i_z).S3 - avS(i_z).T.*r;
    avU(i_z).WrmS3mTr = avU(i_z).W.*r + avU(i_z).S3 + avU(i_z).T.*r;
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
        'Slope',slps(i_v),'Method','direct');
end


% Compensate for plotting

Nlvl = size(avS,1);

for i_z = 1:Nlvl
    r = (1:length(avS(i_z).uuu))*avS(i_z).dr;
    
    for i_v = 1:Nvar
        var = edr_vars{i_v};
        avS(i_z).([var,'_c']) = avS(i_z).(var) ./ r.^slps(i_v) / cons(i_v);
        avU(i_z).([var,'_c']) = avU(i_z).(var) ./ r.^slps(i_v) / cons(i_v);
    end
    
    avS(i_z).edr_S2_4 = avS(i_z).edr_S2*4;
    avU(i_z).edr_S2_4 = avU(i_z).edr_S2*4;
end



%% Save/load

save([myprojectpath,filesep,outputfile],'r_max','transport_method','edr_fit_range','levels','sfc_vars',...
    'S','U','L','N','V','avS','avU','avN',...
    'plotpath')

% load('S_eureca.mat')



%% TABLES

% Dissipation rates

print_vars = {'uu','vv','ww','S2'};%'mS3','WrmS3','WrmS3mTr'};

Nlvl = size(avS,1);
Nvar = numel(print_vars);

fprintf(' %20s',''), fprintf(' & %14s',print_vars{:}), fprintf(' \\\\ \n')
for i_z = 1:Nlvl
    fprintf(' %20s',avS(i_z).level)
    for i_v = 1:Nvar
        var = print_vars{i_v};
        fprintf(' & %6.*f (%.*f)',3,1e4*avS(i_z).(['edr_',var]),3,1e4*avU(i_z).(['edr_',var]) )
    end
    fprintf(' \\\\ \n')
end



%% PLOTS

xlim = [dr r_max];
ylim = [1e-6 5e-3];

Npoints = 40;

Nlvl = size(avS,1);


%% S2

for i_z = 1:Nlvl
    fig = plot_sfc(avS(i_z),avU(i_z),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim,'YLim',[8e-4 1e-2]);
    if i_z>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}C^{-1} $','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
    print(fig,[plotpath,filesep,'s2_',levels{i_z}],'-dpng','-r300')
end
    

%% W

for i_z = 1:Nlvl
    fig = plot_sfc(avS(i_z),avU(i_z),{'W','Wt','Wq'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    if i_z>0
        legend({'$W$','$W_\theta$','$W_q$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
    print(fig,[plotpath,filesep,'W_',levels{i_z}],'-dpng','-r300')
end


%% S3

% for i_z = 1:Nlvl
%     fig = plot_sfc(avS(i_z),avU(i_z),{'S3r'},[1 6],Npoints,'XLim',xlim,'YLim',ylim);
%     if i_z>0
%         legend({'$S_3 r^{-1}$'},'Interpreter','latex','Location','best')
%     end
%     ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
%     print(fig,[plotpath,filesep,'s3r_',levels{i_z}],'-dpng','-r300')
% end


%% T

% for i_z = 1:Nlvl
%     fig = plot_sfc(avS(i_z),avU(i_z),{'T'},5,Npoints,'XLim',xlim,'YLim',ylim);
%     if i_z>0
%         legend({'$T_u$'},'Interpreter','latex','Location','best')
%     end
%     ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
%     print(fig,[plotpath,filesep,'T_',levels{i_z}],'-dpng','-r300')
% end


%% Balance: W, S3/r, T, epsilon

for i_z = 1:Nlvl
    fig = plot_sfc(avS(i_z),avU(i_z),{'W','-S3r','-T','edr_S2_4'},[7 1 5 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_z>0
        legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_z}],'-dpng','-r300')
end

% for i_l = 1:Nlvl
%     fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','-T','edr_S2_4'},[7 1 5 3],Npoints,'XLim',xlim,'YScale','linear');
%     if i_l>0
%         legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
%     end
%     ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
%     print(fig,[plotpath,filesep,'balance_lin_',levels{i_l}],'-dpng','-r300')
% end


%% Total sum: W - S3/r - T

for i_z = 1:Nlvl
    fig = plot_sfc(avS(i_z),avU(i_z),{'WmS3rmT','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_z>0
        legend({'$W-S_3 r^{-1}-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt)) 
    print(fig,[plotpath,filesep,'total_',levels{i_z}],'-dpng','-r300')
end

% for i_l = 1:Nlvl
%     fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3rmT','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YScale','linear');
%     if i_l>0
%         legend({'$W-S_3 r^{-1}-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
%     end
%     ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
%     print(fig,[plotpath,filesep,'total_lin_',levels{i_l}],'-dpng','-r300')
% end


%% S3 components

for i_z = 1:Nlvl
    fig = plot_sfc(avS(i_z),avU(i_z),{'uuu3r','vvu3r','wwu3r'},[7 5 1 2 4 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_z>0
        legend({'uuu','vvu','wwu'},'Interpreter','latex','Location','best')
    end
    ylabel('$3S_3 r^{-1}\,[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
    print(fig,[plotpath,filesep,'s3_comp_',levels{i_z}],'-dpng','-r300')
end


%% Fitting dissipation

for i_z = 1:Nlvl
    fig = plot_sfc_edr(avS(i_z),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim,'YLim',[1e-3 2e-2]);
    if i_z>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_z},avS(i_z).alt))
    print(fig,[plotpath,filesep,'edr2fit_',levels{i_z}],'-dpng','-r300')
end

% for i_l = 1:Nlvl
%     fig = plot_sfc_edr(avS(i_l),[],{'mS3','WrmS3','WrmS3mTr'},[1 7 4],Npoints,'XLim',xlim,'YLim',[1e-5 1e-2]);
%     if i_l>0
%         legend({'$-S_3 r^{-1}$','$W-S_3 r^{-1}$','$W-S_3 r^{-1}-T_u$'},'Interpreter','latex','Location','best')
%     end
%     ylabel('$[\mathrm{m^3\,s^{-3}}]$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
%     print(fig,[plotpath,filesep,'edr3fit_',levels{i_l}],'-dpng','-r300')
% end
