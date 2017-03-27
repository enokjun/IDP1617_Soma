%% Clear Previous File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Creator: Enok Cheon
Date: 10th Mar 2017
Use: IDP deep foundation and male secant pile wall reinforcement design
Key Assumptions:
	- concrete class C30/37
    - steel class S500
    - rectangular stress block

References:
    - Interaction diagrams for reinforced concrete circular cross-section by  & Davor GrandiÄ‡ 

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% factors
gamma_c = 1.5;
gamma_s = 1.15;
alpha_cc = 0.85;

%% Concrete - C30/37 (EC2 TAle 3.1)
% strength 
fck = 30;        % unit: N/mm^2
% dimensions
cover = (25+75);      % unit: mm     % c_min + delta(c_dev)
Col_dia = 800;        % unit: mm
location_pile = 'Riverside Flood Defence - Wall L - non bearing';   % 'edge'

%% Steel Rebar - S500 (EC2 TAle 3.1) 
% strength 
euk = 0.05;
fyk = 500;           % unit: N/mm^2
Es = 200000;        % unit: N/mm^2
% dimensions
Nmin = 3;              % number of longitudinal bars - minimum
Rebar_dia = [16, 20, 25, 32, 40];     % unit: mm     % bar diameter - minimum  
stirrup_dia = [10, 12, 16, 20, 25, 32, 40];   % unit: mm
spacing_max = 200;  % unit: mm     % bar spacing - maximum
spacing_min = 30;   % unit: mm     % bar spacing - minimum - for aggregate sizes

%% loads
% for axial (V): compression (+) and tension (-)
% for moment (M): clockwise (+) and anti-clockwise (-)
% the loads are already factored from the ULS analysis
P_load = -2300;         % unit: kN
M_load = 0;           % unit: kNm
V_load = 0;           % unit: kN

%% analysis
x_increment = 0.1;     % unit: mm

%% plotting
plot_lines = {'b-','g-','r-','k-','m-'}; 
plot_lines_notwork = {'b--','g--','r--','k--','m--'}; 
plot_points = {'bx','gx','rx','kx','mx'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concrete
fcd = (fck*alpha_cc)/gamma_c;     % design concrete yield stress     unit: N/mm^2

% steel
fyd = (fyk/gamma_s);         % design steel yield stress     unit: N/mm^2
esd = fyd/Es;               % design strain of steel

% rectangular stress-strain
ecu = 0.0035;       % max concrete strain
if fck <= 50
    eta = 1;
    lambda = 0.8;
elseif 50 < fck && fck <= 90
    eta = 1-(fck-50)/200;
    lambda = 0.8-(fck-50)/400;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UID - axial and moment - for each diameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_x = cell(1,length(Rebar_dia));
results_y = cell(1,length(Rebar_dia));
N = Nmin;
while N >= Nmin
    d_prime = cover+stirrup_dia(1)+0.5*Rebar_dia(1);
    radius = 0.5*Col_dia - d_prime;
    theta_angle = 2*pi/N;

    spacing_bwt_bars = radius*sin(theta_angle)/sin(0.5*(pi-theta_angle));

    if spacing_bwt_bars < spacing_min || spacing_bwt_bars > spacing_max
        N = N+1;
        continue
    else
        for diameters_to_design = 1:length(Rebar_dia)
            % layout - polar coordinates
            polar_layout = [];
            for num = 1:1:N
                theta_n = 0.5*pi - theta_angle*(num-1);
                polar_layout(num,:) = [radius, theta_n];
            end
            % layout - cartesian coordinates
            [x_coord, y_coord] = pol2cart(polar_layout(:,2), polar_layout(:,1));
            x_coord = round(x_coord,2,'decimal');
            y_coord = round(y_coord,2,'decimal');
            
            % find number of rebars per row/col
            xx = unique(x_coord);
            N_dx = [];
            dx = [];
            for i = 1:length(xx)
                n = find(x_coord == xx(i));
                N_per_d = length(n);
                X_coord = x_coord(n(1));
                N_dx(i) = [N_per_d];
                dx(i) = [X_coord];
            end

            yy = unique(y_coord);
            N_dy = [];
            dy = [];
            for i = 1:length(yy)
                n = find(y_coord == yy(i));
                N_per_d = length(n);
                Y_coord = y_coord(n(1));
                N_dy(i) = [N_per_d];
                dy(i) = [Y_coord];
            end

            % re-configure d -> d1 = furthest and d(last) = closest = d'
            dxx = [];
            dxx(1) = d_prime;
            for j = 2:length(dx)
                dxx(j) = abs(dx(j)-dx(1))+d_prime;
            end
            dxx = flip(dxx);

            dyy = [];
            dyy(1) = d_prime;
            for j = 2:length(dy)
                dyy(j) = abs(dy(j)-dy(1))+d_prime;
            end
            dyy = flip(dyy);

            % area
            Ag = 0.25*pi*Col_dia^2;     % gross area - unit: mm^2
            As_per_bar = 0.25*pi*Rebar_dia(diameters_to_design)^2;     % rebar area - unit: mm^2
            As_dx = As_per_bar*N_dx;
            As_dy = As_per_bar*N_dy;
            As_total = N*As_per_bar;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Longitudinal Reinforcement - check
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reinforcement area check - ec2 9.8.5(3)
            if diameters_to_design < length(Rebar_dia)
                if Ag <= 0.5e6 && As_total < 0.005*Ag
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                elseif (Ag > 0.5e6 && Ag <= 1e6) && As_total < 2500
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                elseif Ag > 1e6 && As_total < 0.0025*Ag
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                else
                    lines = plot_lines{diameters_to_design};
                    points = plot_points{diameters_to_design};
                end
            elseif diameters_to_design == length(Rebar_dia)
                if Ag <= 0.5e6 && As_total < 0.005*Ag
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                elseif (Ag > 0.5e6 && Ag <= 1e6) && As_total < 2500
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                elseif Ag > 1e6 && As_total < 0.0025*Ag
                    lines = plot_lines_notwork{diameters_to_design};
                    points = plot_points{diameters_to_design};
                else
                    lines = plot_lines{diameters_to_design};
                    points = plot_points{diameters_to_design};
                end       
            end
            
            %% axial capacity - compression and tension
            P_comp = (fcd*Ag + As_total*fyd)*1e-3;
            P_ten = -As_total*fyd*1e-3;

            %% X-X iterate by varying the neutral axis (x)
            UIDx = [];
            row = 2;

            % P_comp
            UIDx(1,:) = [0, P_comp];

            %% iterate through different neutral axis
            for x = [Col_dia: -x_increment: 0]
                %% steel - strain, stress and force
                esx = [];
                fsx = [];
                Fsx = [];
                for steel_lines = 1:length(dxx)
                    % strain
                    es_line = ecu*(dxx(steel_lines)-x)/x;

                    % stress 
                    if round(es_line,4,'decimals') == 0
                        fs_line = 0;
                    elseif es_line >= esd 
                        fs_line = fyd;
                    elseif es_line <= -esd
                        fs_line = -fyd;
                    elseif es_line < esd ||  -esd < es_line
                        fs_line = Es*es_line;
                    end

                    % force
                    Fs_line = (As_dx(steel_lines)*fs_line)*1e-3;

                    % tabulate data
                    esx(steel_lines) = es_line;
                    fsx(steel_lines) = fs_line;
                    Fsx(steel_lines) = Fs_line;
                end
                % stop uid when strain is too high for a rebar
                testing_chamber = [];
                for testing = 1:length(esx)
                    if esx(testing) > euk
                        testing_judge = 0;        % if it exceed max steel strain, it fails displayed by 0
                    else 
                        testing_judge = 1;
                    end
                    testing_chamber(testing) = [testing_judge];
                end
                % if all test passed as test_judge with 1, then sum of judge equals the length of matrix
                if length(testing_chamber) ~= sum(testing_chamber)
                    break
                end

                %% concrete - rectangular
                a = lambda*x;
                col_radius = (Col_dia/2);
                if a > col_radius
                    theta = 2*acos((a-col_radius)/col_radius);
                    Ac = Ag - 0.5*(theta-sin(theta))*col_radius^2;
                elseif a == col_radius
                    Ac = Ag/2;
                elseif a < col_radius
                    theta = 2*acos((col_radius-a)/col_radius);
                    Ac = 0.5*(theta-sin(theta))*col_radius^2;
                end
                Fc = (eta*fcd*Ac)*1e-3;

                %% P and M capacity
                P = Fc-sum(Fsx);
                M = (Fc*(col_radius-(a/2)))*1e-3;
                for tension_rebars = 1:length(Fsx)
                    if dxx(tension_rebars) >= col_radius
                        M = M + ((dxx(tension_rebars)-col_radius)*Fsx(tension_rebars))*1e-3;
                    elseif dxx(tension_rebars) < col_radius
                        M = M - ((col_radius-dxx(tension_rebars))*Fsx(tension_rebars))*1e-3;
                    end
                end

                %% tabulate
                UIDx(row,:) = [M, P];

                %% tabulate all data
                row = row + 1;
            end

            % P_ten
            UIDx(row,:) = [0, P_ten];

            %% flter key points
            %{
            key_points_x = [0,P_comp; 0,0; 0,0; 0,P_ten];
            for i = 1:length(UIDx)
                % balance point
                M_rd_max = max(UIDx(:,1));
                if UIDx(i,1) == M_rd_max
                    key_points_x(2,:) = [UIDx(i,1),UIDx(i,2)];
                end
                if  abs(UIDx(i,2)) < 1
                    key_points_x(3,:) = [UIDx(i,1),UIDx(i,2)];
                end
            end
            %}
            %% Plot UID
            figure(1)
            hold on
            plot(UIDx(:,1), UIDx(:,2),lines)
            % plot(key_points_x(:,1),key_points_x(:,2),points)
            grid on

            %% Y-Y iterate by varying the neutral axis (x)
            UIDy = [];
            row = 2;

            % P_comp
            UIDy(1,:) = [0, P_comp];

            %% iterate through different neutral axis
            for x = [Col_dia: -x_increment: 0]
                %% steel - strain, stress and force
                esy = [];
                fsy = [];
                Fsy = [];
                for steel_lines = 1:length(dyy)
                    % strain
                    es_line = ecu*(dyy(steel_lines)-x)/x;

                    % stress 
                    if round(es_line,4,'decimals') == 0
                        fs_line = 0;
                    elseif es_line >= esd 
                        fs_line = fyd;
                    elseif es_line <= -esd
                        fs_line = -fyd;
                    elseif es_line < esd ||  -esd < es_line
                        fs_line = Es*es_line;
                    end

                    % force
                    Fs_line = (As_dy(steel_lines)*fs_line)*1e-3;

                    % tabulate data
                    esy(steel_lines) = es_line;
                    fsy(steel_lines) = fs_line;
                    Fsy(steel_lines) = Fs_line;
                end
                % stop uid when strain is too high for a rebar
                testing_chamber = [];
                for testing = 1:length(esy)
                    if esy(testing) > euk
                        testing_judge = 0;        % if it exceed max steel strain, it fails displayed by 0
                    else 
                        testing_judge = 1;
                    end
                    testing_chamber(testing) = [testing_judge];
                end
                % if all test passed as test_judge with 1, then sum of judge equals the length of matrix
                if length(testing_chamber) ~= sum(testing_chamber)
                    break
                end

                %% concrete - rectangular
                a = lambda*x;
                col_radius = (Col_dia/2);
                if a > col_radius
                    theta = 2*acos((a-col_radius)/col_radius);
                    Ac = Ag - 0.5*(theta-sin(theta))*col_radius^2;
                elseif a == col_radius
                    Ac = Ag/2;
                elseif a < col_radius
                    theta = 2*acos((col_radius-a)/col_radius);
                    Ac = 0.5*(theta-sin(theta))*col_radius^2;
                end
                Fc = (eta*fcd*Ac)*1e-3;

                %% P and M capacity
                P = Fc-sum(Fsy);
                M = (Fc*(col_radius-(a/2)))*1e-3;
                for tension_rebars = 1:length(Fsy)
                    if dyy(tension_rebars) >= col_radius
                        M = M + ((dyy(tension_rebars)-col_radius)*Fsy(tension_rebars))*1e-3;
                    elseif dyy(tension_rebars) < col_radius
                        M = M - ((col_radius-dyy(tension_rebars))*Fsy(tension_rebars))*1e-3;
                    end
                end

                %% tabulate
                UIDy(row,:) = [M, P];

                %% tabulate all data
                row = row + 1;
            end

            %% P_ten
            UIDy(row,:) = [0, P_ten];

            %% flter key points
            %{
            key_points_y = [0,P_comp; 0,0; 0,0; 0,P_ten];
            for i = 1:length(UIDy)
                % balance point
                M_rd_max = max(UIDy(:,1));
                if UIDy(i,1) == M_rd_max
                    key_points_y(2,:) = [UIDy(i,1),UIDy(i,2)];
                end
                if  abs(UIDy(i,2)) < 1
                    key_points_y(3,:) = [UIDy(i,1),UIDy(i,2)];
                end
            end
            %}
            %% Plot UID
            figure(2)
            hold on
            plot(UIDy(:,1), UIDy(:,2),lines)
            % plot(key_points_y(:,1),key_points_y(:,2),points)
            grid on
            
            %% tabulate overall results
            results_x{diameters_to_design} = UIDx;
            results_y{diameters_to_design} = UIDy;
        end            

        % add info the plotted graphs
        figure(1)
        plot(M_load, P_load, 'k*')
        title(strcat('UIDx',{'   '},' diameter =',{' '},int2str(Col_dia),'mm',{'   '},'number of rebars =',{' '},int2str(N)))
        xlabel('Moment [kNm]')
        ylabel('Axial [kN]')
        legend('16mm', '20mm', '25mm','32mm', '40mm', 'loads')
        % legend('16mm', '16mm - key', '20mm', '20mm - key','25mm', '25mm - key', '32mm', '32mm - key', '40mm', '40mm - key','loads')

        figure(2)
        plot(M_load, P_load, 'k*')
        title(strcat('UIDy',{'   '},' diameter =',{' '},int2str(Col_dia),'mm',{'   '},'number of rebars =',{' '},int2str(N)))
        xlabel('Moment [kNm]')
        ylabel('Axial [kN]')
        legend('16mm', '20mm', '25mm','32mm', '40mm', 'loads')
        % legend('16mm', '16mm - key', '20mm', '20mm - key','25mm', '25mm - key', '32mm', '32mm - key', '40mm', '40mm - key','loads')
    end

    % continue to next number of piles 
    contniue_check_pt = input('Are you satisfied and chosen?: (y/n)            ','s');
    if contniue_check_pt == 'y'
        break
    elseif contniue_check_pt == 'n'
        N = N+1;
    else 
        error('please type "y" to show yes or "n" to show no')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shear Force Reinforcement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% choose which orientation and rebar to design
orientation = input('which orientation do you want (1=XX or 2=YY):    ');
rebar_diameter = input('which rebar diameter do you want (16,20,25,32,40):    ');
if rebar_diameter == 16 || rebar_diameter == 20 || rebar_diameter == 25 || rebar_diameter == 32 || rebar_diameter == 40
    As_per_bar = 0.25*pi*rebar_diameter^2;     % rebar area - unit: mm^2
    As_total = N*As_per_bar;
else 
    error('rebar diameter is nto part of the catalog')
end
if orientation == 1
    d = dxx;
    N_d = N_dx;
elseif orientation == 2
    d = dyy;
    N_d = N_dy;
else 
    error('pick 1 or 2')
end
As_d = As_per_bar*N_d;

%% Shear Resistance Reinforcement 
% v1
if fck <= 60
    v1 = 0.6;
elseif fck > 60
    if 0.9-(fck/200) > 0.5
        v1 = 0.9-(fck/200);
    elseif 0.9-(fck/200) < 0.5
        v1 = 0.5;
    end
end
% cot_theta
cot_theta = 2.5; %input('angle of stirrup - 1 ~ 2.5:    ');
% z
z = d(1) - d(length(d));
% V_rd_max
V_rd_max = (Col_dia*v1*fcd*z)/(cot_theta+(1/cot_theta));
if V_load > V_rd_max
    error('try a different diameter or something... V_load > V_rd_max')
end
% As_sw
As_sw = 1000*V_load/(z*fyd*cot_theta);
% sw_max - round down to closest 25mm interval
sw_max = min([20*rebar_diameter, Col_dia, 400]);
sw_max = floor(sw_max/25)*25;

%% minimum shear reinforcement
% rho_w min
rho_w_min = 0.08*sqrt(fck)/fyk;
% Aw_min
Aw_min = rho_w_min*Col_dia*sw_max;

As_sw_design = max(As_sw, (Aw_min/sw_max));

%% design spacing required
shear_reinforcements = [];
row = 1;
for stirrup_diameter = [10, 12, 16, 20, 25, 32, 40]
    for t_spacing = sw_max:-25:50
        As_sw_rd = (2*pi*0.25*stirrup_diameter^2)/t_spacing;
        if As_sw_design > As_sw_rd
            continue
        elseif As_sw_design <= As_sw_rd
            shear_reinforcements(row,:) = [stirrup_diameter, t_spacing, As_sw_rd];
            row=row+1;
        end
    end
end

%% shear reinforcement
[~,order] = sort(shear_reinforcements(:,3));       
shear_reinforcements = shear_reinforcements(order,:);
shear_reinforcements_chosen = shear_reinforcements(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('orientation = ')
disp(orientation)
disp('number of longitudinal bars = ')
disp(N)
disp('longitudinal diameter = ')
disp(rebar_diameter)
disp('As/sw = ')
disp(As_sw_design)
disp('Shear Reinforcement = ')
disp(shear_reinforcements_chosen)

design_summary = [orientation, N, rebar_diameter, As_sw_design, shear_reinforcements_chosen];

dlmwrite(strcat(location_pile,'.csv'),design_summary,'delimiter',',')
disp('done')