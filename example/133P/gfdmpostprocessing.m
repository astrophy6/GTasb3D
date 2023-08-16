% post-processing of GTasb3D code for JGR:Planets MBC study
% written by Yun Zhang

clc;
clear all;

%%% parameter definition
G = 6.67430e-11;
Delta_t = 10;  % timestep, s
m_h2o = 2.9915e-26;  % h2o mass, kg
k_B = 1.380649e-23;  % Boltzmann constant
omega_rot = 2*pi/3600/3.471; % rotation rate, rad/s
alpha_bf = 0.1;  % baskflux fraction
f_d = 0.55; % dust volume fraction
f_i = 0.24; % ice volume fraction
rho_d = 2000;  % dust density
rho_i = 933;   % ice density
rho_full = (f_d*rho_d + f_i*rho_i);
r2i_ratio = f_d*rho_d / (f_i*rho_i);
flag_gas = 1;  % 0: for simulations without gas; 1: for simulation with gas
r_d = 0.001; % average dust particle radius, m
F_cohesion = 0.1e-6;  % cohesive force, N
D_max = 1;  % upper limit of the SFD, m
Volume_D_max = Volume_SFD(D_max);
N_ice_layer = 9; % the deepest possible layer of the ice-sublimation layer

%%% input parameter control - change as needed
nInitial = 17760000; 
nOutInterval = 80000; 
ntotal = 35520000;
filedir = "dustmanle0cm_obliquity75_ExampleOutcomes"; % change to the directory name as needed

%%% read node data
fileID = fopen('133P_ellipsoid.txt','r');
num_node = floor(fread(fileID,1,'double'));
nLayer1D = floor(fread(fileID,1,'double'));
dArea_1Dlayer = fread(fileID,[1,nLayer1D],'double'); 
node_data = fread(fileID,[4,num_node],'double'); 
node_data = node_data';
fclose(fileID);
N_layer = sum(node_data(:,1)==0);
for i = 0:20
    depth(i+1) = norm(node_data(N_layer*i+1,2:4)-node_data(N_layer*(i+1)+1,2:4));
    if depth(i+1) > 100.0
        disp(['1D node depth: ', num2str(sum(depth(1:(i-1)))), ' m']);
        depth(i+1) = norm(node_data(N_layer*i+1,2:4)) - max(vecnorm(node_data((N_layer*(i+1)+1):(N_layer*(i+2)+1),2:4),2,2));
        break
    end
end
A_sur = dArea_1Dlayer(1);  % total surface area

%%% prepare gravity calculation based on shape information
axis_a = 2700;
axis_b = 1800;
axis_c = 1800;
alpha = axis_c/axis_a;
beta = axis_b/axis_a;
syms fun_g u_g
fun_g = 1.0 / ( (u_g+1.0)^(3/2) * (u_g+beta^2)^(1/2) * (u_g+alpha^2)^(1/2) );
A_x_syms = alpha * beta * vpaintegral(fun_g,u_g,0,inf);
A_x = double(A_x_syms);
fun_g = 1.0 / ( (u_g+1.0)^(1/2) * (u_g+beta^2)^(3/2) * (u_g+alpha^2)^(1/2) );
A_y_syms = alpha * beta * vpaintegral(fun_g,u_g,0,inf);
A_y = double(A_y_syms);
fun_g = 1.0 / ( (u_g+1.0)^(1/2) * (u_g+beta^2)^(1/2) * (u_g+alpha^2)^(3/2) );
A_z_syms = alpha * beta * vpaintegral(fun_g,u_g,0,inf);
A_z = double(A_z_syms);
f_body(1,:) = ( - 2*pi*rho_full*G*A_x + omega_rot^2 ) * node_data(:,2)';
f_body(2,:) = ( - 2*pi*rho_full*G*A_y + omega_rot^2 ) * node_data(:,3)';
f_body(3,:) = ( - 2*pi*rho_full*G*A_z ) * node_data(:,4)';


%%% read thermophysical variable data and calculate the gas and dust production rate
i = 1;
for num_frame = nInitial:nOutInterval:ntotal
    filename = filedir+"/gfdm."+sprintf('%09d',num_frame);
    if flag_gas
        [time,f,data,data_gas] = ReadGfdm(filename,num_node,flag_gas,node_data);
    else
        [time,f,data] = ReadGfdm(filename,num_node,flag_gas,node_data);
    end
    Time(i) = time;
    f_all(i) = f;
    temp(i,:) = data(:,1)';
    Rho_gas(i,:) = data(:,2)';
    f_ice(i,:) = data(:,3)';
    p_gas(i,:) = Rho_gas(i,:) .* temp(i,:) * k_B/m_h2o;
    mass_rate(i,:) = (1-alpha_bf)*Rho_gas(i,node_data(:,1)==0).*sqrt(temp(i,node_data(:,1)==0)*k_B/2/pi/m_h2o);
    Q_rate_all(i) = sum(mass_rate(i,:))*A_sur/N_layer;

    % find ice-sublimation front
    for j = 1:N_layer
        for k = 1:N_ice_layer
            if f_ice(i, k*N_layer+j) > 1e-6  % ice-containing
                N_front(i,j) = k*N_layer+j;
                N_front_layer(i,j) = k;
                break;
            end
        end
    end

    % dust levitation
    f_drag = pi*p_gas(i,:);
    for j = 1:N_layer
        norm_surface_0 = node_data(j+N_layer,2:4) - node_data(j,2:4);
        norm_surface = norm_surface_0/norm(norm_surface_0);
        f_body_norm(j) = norm_surface*f_body(:,N_front(i,j)) * 4/3*pi*rho_d;
        f_hydro_pres = 0;
        poly_fun = [f_body_norm(j), f_hydro_pres-f_drag(N_front(i,j)), 0, F_cohesion];
        r_root = roots(poly_fun);
        j_real = 1;
        for jj = 1:3
            if isreal(r_root(jj)) && r_root(jj) > 0
                r_root_real(j_real) = r_root(jj);
                j_real = j_real + 1;
            end
        end
        if j_real > 1
            r_dust_min(i,j) = min(r_root_real);
            r_dust_max_temp = max(r_root_real);
            r_dust_max(i,j) = min(r_dust_max_temp, 0.5*sum(depth(1:N_front_layer(i,j))));

            if r_dust_max(i,j) > r_dust_min(i,j)
                f_flux_d(i,j) = mass_rate(i,j) * r2i_ratio * ( Volume_SFD(2*r_dust_max(i,j)) - Volume_SFD(2*r_dust_min(i,j)) ) / Volume_D_max;
            else
                r_dust_min(i,j) = 0.0;
                r_dust_max(i,j) = 0.0;
                f_flux_d(i,j) = 0.0;
            end
        else
            r_dust_min(i,j) = 0.0;
            r_dust_max(i,j) = 0.0;
            f_flux_d(i,j) = 0.0;
        end
    end
    i = i + 1;
end
i_last = i - 1;


%%% max surface evolution 
[~,theta_sur,~] = cart2sph(node_data(node_data(:,1)==0,2),node_data(node_data(:,1)==0,3),node_data(node_data(:,1)==0,4));
for i = 1:size(temp,1)
    for j = 1:73
        lat_inter = (j-37)/72*pi;
        [temp_max(i,j),temp_max_index(i,j)] = max(temp(i,abs(theta_sur-lat_inter)<pi/72));
        gas_rate_max(i,j) = max(mass_rate(i,abs(theta_sur-lat_inter)<pi/72));
        dust_r_max(i,j) = max(r_dust_max(i,abs(theta_sur-lat_inter)<pi/72));
        dust_r_min(i,j) = max(r_dust_min(i,abs(theta_sur-lat_inter)<pi/72));
        dust_rate_max(i,j) = max(f_flux_d(i,abs(theta_sur-lat_inter)<pi/72));
    end
end
[XX,YY] = meshgrid(f_all',-90:2.5:90);

figure;
surf((XX-3*pi)/pi*180+180,YY,temp_max','EdgeAlpha',0.0);
box on;
shading interp;
colormap pink; hc = colorbar; hc.Label.String = 'Diurnal maximum temperature $T$ [K]';
ylabel('Latitude');
xlabel('True anomaly');
xlim([0,360]+180);
ylim([-90,90]);
set(gcf,'Position',[10 10 650 170]);
set(gca,'Position',[0.13 0.31 0.67 0.62],'FontName','Times New Roman','FontSize',22);
view(0,90);
xticks([0, 90, 180, 270, 360]+180);
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-90^\circ','-45^\circ','0^\circ','45^\circ','90^\circ'});
yticks([-90,-45,0,45,90]);
caxis([120,240]);

figure;
surf((XX-3*pi)/pi*180+180,YY,log10(gas_rate_max'),'EdgeAlpha',0.0);
box on;
shading interp;
colormap jet; hc = colorbar; hc.Label.String = 'Diurnal maximum gas number density log10($n_\mathrm{g}$) [$\mathrm{m}^{-3}$]';
ylabel('Latitude');
xlabel('True anomaly');
xlim([0,360]+180);
ylim([-90,90]);
set(gcf,'Position',[10 10 650 170]);
set(gca,'Position',[0.13 0.31 0.67 0.62],'FontName','Times New Roman','FontSize',22);
view(0,90);
xticks([0, 90, 180, 270, 360]+180);
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-90^\circ','-45^\circ','0^\circ','45^\circ','90^\circ'});
yticks([-90,-45,0,45,90]);
caxis([-10,-5]);

figure;
surf((XX-3*pi)/pi*180+180,YY,log10(dust_r_max'/1e-3),'EdgeAlpha',0.0);
box on;
shading interp;
colormap jet; hc = colorbar; hc.Label.String = 'Diurnal largest liftable dust particle radius log10($r_\mathrm{d}$) [mm]';
ylabel('Latitude');
xlabel('True anomaly  [deg]');
xlim([0,360]+180);
ylim([-90,90]);
set(gcf,'Position',[10 10 650 170]);
set(gca,'Position',[0.13 0.31 0.67 0.62],'FontName','Times New Roman','FontSize',22);
view(0,90);
xticks([0, 90, 180, 270, 360]+180);
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-90^\circ','-45^\circ','0^\circ','45^\circ','90^\circ'});
yticks([-90,-45,0,45,90]);
caxis(log10([0.05,500]));

figure;
surf((XX-3*pi)/pi*180+180,YY,log10(dust_r_min'/1e-3),'EdgeAlpha',0.0);
box on;
shading interp;
colormap jet; hc = colorbar; hc.Label.String = 'Diurnal largest liftable dust particle radius log10($r_\mathrm{d}$) [mm]';
ylabel('Latitude');
xlabel('True anomaly  [deg]');
xlim([0,360]+180);
ylim([-90,90]);
set(gcf,'Position',[10 10 650 170]);
set(gca,'Position',[0.13 0.31 0.67 0.62],'FontName','Times New Roman','FontSize',22);
view(0,90);
xticks([0, 90, 180, 270, 360]+180);
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-90^\circ','-45^\circ','0^\circ','45^\circ','90^\circ'});
yticks([-90,-45,0,45,90]);
caxis(log10([0.05,500]));

figure;
surf((XX-3*pi)/pi*180+180,YY,log10(dust_rate_max)','EdgeAlpha',0.0);
box on;
shading interp;
colormap parula; hc = colorbar; hc.Label.String = 'Diurnal maximum dust mass loss flux log10($\dot{m}_\mathrm{d}$) [$\mathrm{kg~m}^{-2}~\mathrm{s}^{-1}$]';
ylabel('Latitude');
xlabel('True anomaly  [deg]');
xlim([0,360]+180);
ylim([-90,90]);
set(gcf,'Position',[10 10 650 170]);
set(gca,'Position',[0.13 0.31 0.67 0.62],'FontName','Times New Roman','FontSize',22);
view(0,90);
xticks([0, 90, 180, 270, 360]+180);
xticklabels({'-180','-90','0','90','180'});
yticklabels({'-90^\circ','-45^\circ','0^\circ','45^\circ','90^\circ'});
yticks([-90,-45,0,45,90]);
caxis([-8,-5]);


%% production rate verse observation
dust_rate_obs = 1.4; % kg/s (Jewitt et al. 2014)
Y_data = [-90:2.5:90];
[~,theta_sur,~] = cart2sph(node_data(node_data(:,1)==0,2),node_data(node_data(:,1)==0,3),node_data(node_data(:,1)==0,4));
Area_element = A_sur/N_layer;

%%% impact
lat_crater_index = [1,7,13,19,25,31,37,43,49,55,61,67,73]; % [-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90]
crater_index_longitude0 = [3771, 1397, 790, 1391, 3704, 2225, 941, 2683, 1272, 278, 3501, 3229, 337];  % crater center surface node No.
crater_index_longitude90 = [3771, 3341, 3622, 3697, 1574, 3353, 2048, 3639, 3224, 3060, 698, 3441, 337];

for i = 1:length(lat_crater_index)
    crater_center0 = node_data(crater_index_longitude0(i),2:4);
    crater_center90 = node_data(crater_index_longitude90(i),2:4);
    crater_elements0 = [crater_index_longitude0(i)];
    crater_elements90 = [crater_index_longitude90(i)];
    Q_dust_crater0(i,:) = sum(f_flux_d(:,crater_elements0)',1)*Area_element;
    Q_dust_crater90(i,:) = sum(f_flux_d(:,crater_elements90)',1)*Area_element;

    node_distance0 = vecnorm(node_data(1:N_layer,2:4)-crater_center0, 2, 2);
    [distance0, index0] = sort(node_distance0);
    node_distance90 = vecnorm(node_data(1:N_layer,2:4)-crater_center90, 2, 2);
    [distance90, index90] = sort(node_distance90);
    next_closest0 = 2;
    next_closest90 = 2;

    while max(Q_dust_crater0(i,:)) < dust_rate_obs
        crater_elements0 = [crater_elements0, index0(next_closest0)];
        next_closest0 = next_closest0 + 1;
        Q_dust_crater0(i,:) = sum(f_flux_d(:,crater_elements0)',1)*Area_element;
    end

    while max(Q_dust_crater90(i,:)) < dust_rate_obs
        crater_elements90 = [crater_elements90, index90(next_closest90)];
        next_closest90 = next_closest90 + 1;
        Q_dust_crater90(i,:) = sum(f_flux_d(:,crater_elements90)',1)*Area_element;
    end

    r_crater0(i) = sqrt(size(crater_elements0,2)*Area_element/pi);
    r_crater90(i) = sqrt(size(crater_elements90,2)*Area_element/pi);
    Q_gas_crater0(i,:) = sum(mass_rate(:,crater_elements0)',1)*Area_element;
    Q_gas_crater90(i,:) = sum(mass_rate(:,crater_elements90)',1)*Area_element;
    Q_gas_total0(i) = sum(Q_gas_crater0(i,:))*Delta_t*nOutInterval;
    Q_dust_total0(i) = sum(Q_dust_crater0(i,:) )*Delta_t*nOutInterval;
    Q_gas_total90(i) = sum(Q_gas_crater90(i,:) )*Delta_t*nOutInterval;
    Q_dust_total90(i) = sum(Q_dust_crater90(i,:) )*Delta_t*nOutInterval;
end

figure;
subplot(2,2,1);
hold on;
box on;
for i = 1:length(lat_crater_index)
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater0(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_dust_crater0(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_dust_crater0(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total dust production rate [kg/s]');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.55 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'','','','',''});
yticks([0,0.4,0.8,1.2,1.4,1.6]);
subplot(2,2,3);
hold on;
box on;
for i = 1:length(lat_crater_index)
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater0(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_gas_crater0(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_gas_crater0(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total gas production rate [kg/s]');
xlabel('True anomaly');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.12 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'-180','-90','0','90','180'});
yticks([0,0.4,0.8,1.2,1.6]);

subplot(2,2,2);
hold on;
box on;
for i = 1:length(lat_crater_index)
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater90(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_dust_crater90(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_dust_crater90(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total dust production rate [kg/s]');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.6 0.55 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'','','','',''});
yticks([0,0.4,0.8,1.2,1.4,1.6]);
subplot(2,2,4);
hold on;
box on;
for i = 1:length(lat_crater_index)
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater90(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_gas_crater90(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_gas_crater90(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total gas production rate [kg/s]');
xlabel('True anomaly');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.6 0.12 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'-180','-90','0','90','180'});
yticks([0,0.4,0.8,1.2,1.6]);



%%% impact with 249-m-radius crater
Num_element = 14;
for i = 1:length(lat_crater_index)
    crater_center0 = node_data(crater_index_longitude0(i),2:4);
    crater_elements0 = [crater_index_longitude0(i)];
    Q_dust_crater0(i,:) = sum(f_flux_d(:,crater_elements0)',1)*Area_element;

    node_distance0 = vecnorm(node_data(1:N_layer,2:4)-crater_center0, 2, 2);
    [distance0, index0] = sort(node_distance0);
    next_closest0 = 2;

    Num_active = 1;
    while Num_active < Num_element
        crater_elements0 = [crater_elements0, index0(next_closest0)];
        next_closest0 = next_closest0 + 1;
        Q_dust_crater0(i,:) = sum(f_flux_d(:,crater_elements0)',1)*Area_element;
        Num_active = Num_active + 1;
    end

    r_crater0(i) = sqrt(size(crater_elements0,2)*Area_element/pi);
    Q_gas_crater0(i,:) = sum(mass_rate(:,crater_elements0)',1)*Area_element;
    Q_gas_total0(i) = sum(Q_gas_crater0(i,:))*Delta_t*nOutInterval;
    Q_dust_total0(i) = sum(Q_dust_crater0(i,:) )*Delta_t*nOutInterval;
end

figure;
subplot(2,2,1);
hold on;
box on;
for i = 1:length(lat_crater_index)
    if r_crater0(i) == Inf
        continue;
    end
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater0(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_dust_crater0(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_dust_crater0(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total dust production rate [kg/s]');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 300 400]);
set(gca,'Position',[0.12 0.55 0.89 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'','','','',''});
yticks([0,0.4,0.8,1.2,1.4,1.6]);
subplot(2,2,3);
hold on;
box on;
for i = 1:length(lat_crater_index)
    if r_crater0(i) == Inf
        continue;
    end
    str = '$'+string(Y_data(lat_crater_index(i)))+'^\circ$ ($R_\mathrm{crater} = '+ string(round(r_crater0(i))) +'$ m)';
    if Y_data(lat_crater_index(i)) > 0
        plot((f_all-3*pi)/pi*180, Q_gas_crater0(i,:),'LineWidth',2.0, 'DisplayName',str,'linestyle','--');
    else
        plot((f_all-3*pi)/pi*180, Q_gas_crater0(i,:),'LineWidth',2.0, 'DisplayName',str);
    end
end
ylabel('Total gas production rate [kg/s]');
xlabel('True anomaly');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 300 400]);
set(gca,'Position',[0.12 0.12 0.89 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'-180','-90','0','90','180'});
yticks([0,0.4,0.8,1.2,1.6]);


%%% landslides
landslides_percentage = 0.0255;
load('active_surface_255_sol0.mat');
i = 1;
active_percentage = sum(active_surface_255)/N_layer;
Q_gas_landslide(i,:) = sum(mass_rate(:,active_surface_255)')*A_sur/N_layer;
Q_dust_landslide(i,:) = sum(f_flux_d(:,active_surface_255)')*A_sur/N_layer;
Q_gas_landslide_total(i) = sum(Q_gas_landslide(i,:) )*Delta_t*nOutInterval;
Q_dust_landslide_total(i) = sum(Q_dust_landslide(i,:) )*Delta_t*nOutInterval;

sqrt(active_percentage*A_sur/pi)

figure;
subplot(2,1,1);
hold on;
for i = 1:length(landslides_percentage)
    str = string(active_percentage(i)*100)+'% resurfacing';
    plot((f_all-3*pi)/pi*180, Q_dust_landslide(i,:),'LineWidth',2.0, 'DisplayName',str);
end
ylabel('Total dust production rate [kg/s]');
grid on;
xlim([0,360]);
ylim([0,1.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.55 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'','','','',''});
yticks([0,0.4,0.8,1.2,1.4,1.8]);
box on;

subplot(2,1,2);
hold on;
for i = 1:length(landslides_percentage)
    str = string(active_percentage(i)*100)+'% resurfacing';
    plot((f_all-3*pi)/pi*180, Q_gas_landslide(i,:),'LineWidth',2.0, 'DisplayName',str);
end
ylabel('Total gas production rate [kg/s]');
xlabel('True anomaly');
grid on;
xlim([0,360]);
ylim([0,0.6]);
set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.12 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'-180','-90','0','90','180'});
yticks([0,0.4,0.8,1.2,1.4,1.8]);
box on;

% active range: [-10,109] deg
plot([170,170],[0,2.0],'LineWidth',1.0,'Color','r');
plot([289,289],[0,2.0],'LineWidth',1.0,'Color','r');


%%% entire surface exposure
Q_gas_whole = sum(mass_rate')*Area_element;
Q_dust_whole = sum(f_flux_d')*Area_element;
Q_gas_whole_total = sum(Q_gas_whole(i,:) )*Delta_t*nOutInterval;
Q_dust_whole_total = sum(Q_dust_whole(i,:) )*Delta_t*nOutInterval;

figure;
subplot(2,1,1);
hold on;
str = '$\varepsilon_\mathrm{O}=0^\circ$';
plot((f_all-3*pi)/pi*180, Q_dust_whole,'LineWidth',2.0, 'DisplayName',str);
ylabel('Total dust production rate [kg/s]');
grid on;
xlim([0,360]);

set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.55 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'','','','',''});
box on;

subplot(2,1,2);
hold on;
plot((f_all-3*pi)/pi*180, Q_gas_whole,'LineWidth',2.0, 'DisplayName',str);
ylabel('Total gas production rate [kg/s]');
xlabel('True anomaly');
grid on;
xlim([0,360]);

set(gcf,'Position',[10 10 600 400]);
set(gca,'Position',[0.12 0.12 0.39 0.39],'FontName','Times New Roman','FontSize',18);
view(0,90);
xticks([0, 90, 180, 270, 360]);
xticklabels({'-180','-90','0','90','180'});
box on;
