
datapath = '/mnt/mdisk/open/BOMEX/simulation_202407/';

addpath(genpath(myprojectpath))

plotpath = [myprojectpath,filesep,'figures_bomex'];
if ~isfolder(plotpath), mkdir(plotpath), end


vars = {'thl','ql','qt','u','v','w','p'};

t_start = 3*3600+1;
dx_q = 100;
dy_q = 100;


% Load coordinates

t = ncread([datapath,'thl.nc'],'time');
x = ncread([datapath,'thl.nc'],'x');
y = ncread([datapath,'thl.nc'],'y');
z = ncread([datapath,'thl.nc'],'z');
xh = ncread([datapath,'u.nc'],'xh');
yh = ncread([datapath,'v.nc'],'yh');
zh = ncread([datapath,'w.nc'],'zh');
dx = unique(diff(x));
dy = unique(diff(y));
dz = unique(diff(z));


% Load layers

prof = load([myprojectpath,filesep,'bomex_profiles.mat']);


% Define slices

ind_t = find(t>=t_start);
t_q = t(ind_t);

x_q = dx_q:dx_q:x(end);
ind_x = unique(interp1(x,1:length(x),x_q,'nearest'));
x_q = x(ind_x);

y_q = dy_q:dy_q:y(end);
ind_y = unique(interp1(y,1:length(y),y_q,'nearest'));
y_q = y(ind_y);

levels = {'near-surface','mid-subcloud','top-subcloud','cloud-base'};
z_q = [60*ones(size(prof.cb_cf)); 2/3*prof.ml_qv; 0.5*(prof.ml_qv+prof.cb_cf); prof.cb_cf];
z_q = z_q(:,ind_t);
ind_z = interp1(prof.z,1:numel(prof.z),z_q,'nearest');
z_q = z(ind_z);



Nx = length(x_q);
Ny = length(y_q);
[Nz,Nt] = size(ind_z);
Nvar = numel(vars);

TURB = cell(Nz,Nt);
errors = struct([]);

for i_t = 1:Nt
    fprintf('i_t = %2d of %2d: i_z = ',i_t,Nt)
    
    for i_z = 1:Nz
        fprintf(' %d',i_z)
        
        % Load variables in a horizontal plane
        
        plane = cell2struct(cell(Nvar,1),vars);
        
        for i_v = 1:Nvar
            var = vars{i_v};
            if strcmp(var,'w'), zn = 2;
            else, zn = 1; end
            try
                plane.(var) = ncread([datapath,var,'.nc'],var,[1 1 ind_z(i_z,i_t) ind_t(i_t)],[Inf Inf zn 1]);
            catch
                errors = vertcat(errors, struct('i_t',i_t,'i_z',i_z,'var',var));
                plane.(var) = nan(length(x),length(y),zn);
            end
        end
    
        
        % Interpolate to basic grid
        
        plane.u = interp1( [xh;xh(end)+dx], [plane.u;plane.u(1,:)] ,x);
        plane.v = shiftdim( interp1( [yh;yh(end)+dy], shiftdim([plane.v,plane.v(:,1)],1), y) ,1);
        plane.w = shiftdim( interp1( zh(ind_z(i_z,i_t)+[0 1]), shiftdim(plane.w,2), z(ind_z(i_z,i_t))) ,1);
        
        
        % Cut out virtual flight segments
        
        turbx = struct([]);
        for i_x = 1:Nx
            turbx(i_x).level = string(levels{i_z});    
            turbx(i_x).dir   = "cross";
            turbx(i_x).ind_x = ind_x(i_x);      
            turbx(i_x).ind_y = nan;
            turbx(i_x).ind_z = ind_z(i_z,i_t);
            for i_v = 1:Nvar
                var = vars{i_v};
                turbx(i_x).(var) = plane.(var)(ind_x(i_x),:)';
            end 
            turbx(i_x).pbase = prof.p(ind_z(i_z),ind_t(i_t));
            turbx(i_x).rho = prof.rho(ind_z(i_z),ind_t(i_t));
        end
        
        turby = struct([]);
        for i_y = 1:Ny
            turby(i_y).level = string(levels{i_z});    
            turby(i_y).dir   = "along";
            turby(i_y).ind_x = nan;
            turby(i_y).ind_y = ind_y(i_y);    
            turby(i_y).ind_z = ind_z(i_z,i_t);
            for i_v = 1:Nvar
                var = vars{i_v};
                turby(i_y).(var) = plane.(var)(:,ind_y(i_y));
            end
            turby(i_y).pbase = prof.p(ind_z(i_z),ind_t(i_t));
            turby(i_y).rho = prof.rho(ind_z(i_z),ind_t(i_t));
        end
        
        TURB{i_z,i_t} = [turbx';turby'];
        
        
        % Plot plane
        
%         var = 'w';        
%         fig = figure('Color','white','PaperUnits','centimeters',...
%             'PaperSize',[16 9],'PaperPosition',[0 0 16 9]);
%         ax = axes('Box','off','Color','none','FontSize',10,'YDir','normal');
%         hold on
%         imagesc(x,y,plane.(var)')
%         axis tight
%         colorbar
% %         for i_x = 1:Nx
% %             plot(x_q(i_x)*[1 1],ax.YLim,'r','LineWidth',1)
% %         end
% %         for i_y = 1:Nx
% %             plot(ax.XLim,y_q(i_y)*[1 1],'r','LineWidth',1)
% %         end
%         title(sprintf('%s at z=%.0fm t=%.2fh',var,z_q(i_z),t_q(i_t)/3600))
%         print(fig,[plotpath,filesep,sprintf('plane_%s_z%04d_t%03d',var,z_q(i_z),t_q(i_t)/60)],'-r300','-dpng')
        
    end
    
    fprintf('\n')
    close all
    
end

for i_e = 1:length(errors)
    fprintf('Unable to load i_t=%d, i_z=%d, t=%.0f, z=%.0f, var=%s\n',...
        errors(i_e).i_t,errors(i_e).i_z,t_q(errors(i_e).i_t),z_q(errors(i_e).i_z),errors(i_e).var)
    TURB{errors(i_e).i_z,errors(i_e).i_t} = [];
end

TURB = vertcat(TURB{:});

save('bomex_turb.mat','TURB','t','x','y','z','dx','dy','dz')

