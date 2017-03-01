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
Date: 18th Feb 2017
Use: IDP embedded wall reinforcement design
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
fctm = 2.9;		 % unit: N/mm^2
% dimensions
cover = 25+10;     	 % unit: mm     % c_min + delta(c_dev)
depths = 600;        % unit: mm
widths = 600;		 % unit: mm

%% Steel Rebar - S500 (EC2 TAle 3.1) 
% strength 
euk = 0.05;
fyk = 500;           % unit: N/mm^2
Es = 200000;        % unit: N/mm^2
% dimensions
N_max = 6;              % number of longitudinal bars - flexure - max
N_max_tor = 8;			% number of longitudinal bars - torsion - max
Rebar_dia = [6, 8, 10, 12, 16, 20, 25, 32, 40];     % unit: mm     % bar diameter - minimum  
stirrup_diameter = [8, 10, 12, 16, 20, 25, 32, 40];   % unit: mm
spacing_max = 200;  % unit: mm     % bar spacing - maximum
spacing_min = 30;   % unit: mm     % bar spacing - minimum - for aggregate sizes

%% loads
% for axial (V): compression (+) and tension (-)
% for moment (M): clockwise (+) and anti-clockwise (-)
% the loads are already factored from the ULS analysis
M_pos_load = 120;           % Moment+   unit: kNm
M_neg_load = 100;           % Moment-   unit: kNm
V_load = 100;               % Shear     unit: kN
T_load = 50;                % Torsion   unit: kNm

%% analysis
x_increment = 0.1;     % unit: mm
K0 = 0.168;

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
%% Bending Moment Design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Positive Bending Moment 
d_pos = depths - cover - (stirrup_diameter(1)/2);
d_prime_pos = cover+(stirrup_diameter(1)/2);
K_pos = (1e06)*M_pos_load/(fck*widths*d_pos^2);
if K_pos <= K0
	Z = d_pos*(0.5+sqrt(0.25-(3/3.4)*K_pos));
	x = (Z-d_pos)/0.4;
	if x > 0.45*d_pos
		error('nuetal axis is too low')
	end
	As_pos_1 = (1e06)*M_pos_load/(0.87*fyk*Z);
	As_neg_1 = 0;
elseif K_pos > K0
	M_rd_max = (1e-06)*K0*fck*widths*d_pos^2;
	if d_prime_pos/(0.45*d_pos) > 0.38
		error('decrease d_prime_pos')
	end
	As_neg_1 = (1e06)*(M_pos_load-M_rd_max)/(0.87*fyk*(d_pos-d_prime_pos));
	Z = d_pos*(0.5+sqrt(0.25-(3/3.4)*K0));
	As_pos_1 = (1e06)*(M_pos_load/(0.87*fyk*Z))+As_neg_1;
end

%% Negative Bending Moment 
d_neg = depths - cover - (stirrup_diameter(1)/2);
d_prime_neg = cover+(stirrup_diameter(1)/2);
K_neg = (1e06)*M_neg_load/(fck*widths*d_neg^2);
if K_neg <= K0
	Z = d_neg*(0.5+sqrt(0.25-(3/3.4)*K_neg));
	x = (Z-d_neg)/0.4;
	if x > 0.45*d_neg
		error('nuetal axis is too low')
	end
	As_neg_2 = (1e06)*M_pos_load/(0.87*fyk*Z);
	As_pos_2 = 0;
elseif K_neg > K0
	M_rd_max = (1e-06)*K0*fck*widths*d_neg^2;
	if d_prime_neg/(0.45*d_pos) > 0.38
		error('decrease d_prime_neg')
	end
	As_pos_2 = (1e06)*(M_neg_load-M_rd_max)/(0.87*fyk*(d_neg-d_prime_neg));
	Z = d_pos*(0.5+sqrt(0.25-(3/3.4)*K0));
	As_neg_2 = (1e06)*(M_neg_load/(0.87*fyk*Z))+As_pos_2;
end

%% Required Postive and negative reinforcement areas
As_pos = max(As_pos_1, As_pos_2);
As_neg = max(As_neg_1, As_neg_2);

rebar_areas_pos = [];
rebar_areas_neg = [];
row_p = 1;
row_n = 1;
for diameter_rebar = Rebar_dia
	As_per_bar = 0.25*pi*diameter_rebar^2;
	n_pos = ceil(As_pos/As_per_bar);
	n_neg = ceil(As_neg/As_per_bar);
	A_pos_rd = As_per_bar*n_pos;
	A_neg_rd = As_per_bar*n_neg;
    if n_pos <= N_max
        rebar_areas_pos(row_p,:) = [As_pos, diameter_rebar, As_per_bar, n_pos, A_pos_rd];
        row_p = row_p+1;
    end
    if n_neg <= N_max
        rebar_areas_neg(row_n,:) = [As_neg, diameter_rebar, As_per_bar, n_neg, A_neg_rd];
        row_n = row_n+1;
    end
end

[~,order] = sort(rebar_areas_pos(:,5));       
rebar_areas_pos = rebar_areas_pos(order,:);
[~,order] = sort(rebar_areas_neg(:,5));       
rebar_areas_neg = rebar_areas_neg(order,:);

% final reinforcement area
rebar_areas_pos_chosen = rebar_areas_pos(1,:);
rebar_areas_neg_chosen = rebar_areas_neg(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check - longitudinal area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = depths-cover-(stirrup_diameter(1)/2);
As_min = max(0.26*(fctm/fyk)*widths*d, 0.0013*widths*d);
As_max = 0.04*depths*widths;

if rebar_areas_pos_chosen(5) < As_min 
	new_As = find(rebar_areas_pos(:,5) >= As_min);
	rebar_areas_pos_chosen = rebar_areas_pos(new_As(1),:);
elseif rebar_areas_pos_chosen(5) > As_max
	new_As = find(rebar_areas_pos(:,5) <= As_max);
	rebar_areas_pos_chosen = rebar_areas_pos(new_As(1),:);
end

if rebar_areas_neg_chosen(5) < As_min 
	new_As = find(rebar_areas_neg(:,5) >= As_min);
	rebar_areas_neg_chosen = rebar_areas_neg(new_As(1),:);
elseif rebar_areas_neg_chosen(5) > As_max
	new_As = find(rebar_areas_neg(:,5) <= As_max);
	rebar_areas_neg_chosen = rebar_areas_neg(new_As(1),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shear Force Reinforcement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% concrete shear design
d = depths-cover-(stirrup_diameter(1)/2);
Crd_c = 0.18/gamma_c;	
k = min(1+sqrt(200/d),2);		
rho1 = min(rebar_areas_pos(5)/(widths*d),0.02);	
vmin = 0.035*(k^1.5)*sqrt(fck);
Vrd_c = max(Crd_c*k*((100*rho1*fck)^(1/3))*widths*d,vmin*widths*d);

%% Shear Resistance Reinforcement/Torsional Reinforcement 
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
z = 0.9*d;
% V_rd_max
V_rd_max = (widths*v1*fcd*z)/(cot_theta+(1/cot_theta));
if V_load > V_rd_max
    error('try a different diameter or something... V_load > V_rd_max')
end

% Does it require steel reinforcement?
if Vrd_c > V_load
	%% minimum shear reinforcement
	% sw_max
	sw_max = min([20*min(rebar_areas_pos(2),rebar_areas_neg(2)), min(widths,depths), 400]);
	% rho_w min
	rho_w_min = 0.08*sqrt(fck)/fyk;
	% Aw_min
	Aw_min = rho_w_min*widths*sw_max;
	As_sw = 0;
	
	%% design spacing required
	shear_reinforcements = [];
	row = 1;
	for stirrup_dia = stirrup_diameter
        As_rd = (2*pi*0.25*stirrup_dia^2);
        As_sw_rd = As_rd/sw_max;
        if As_rd < Aw_min
            continue
        elseif As_rd >= Aw_min
            shear_reinforcements(row,:) = [stirrup_dia, sw_max, As_sw_rd];
            row=row+1;
        end
	end
elseif Vrd_c <= V_load
	% As_sw
	As_sw = 1000*V_load/(z*fyd*cot_theta);
	% sw_max
	sw_max = min([20*min(rebar_areas_pos(2),rebar_areas_neg(2)), min(widths,depths), 400]);

	%% design spacing required
	shear_reinforcements = [];
	row = 1;
	for stirrup_dia = stirrup_diameter
	    for t_spacing = sw_max:-25:50
	        As_sw_rd = (2*pi*0.25*stirrup_dia^2)/t_spacing;
	        if As_sw > As_sw_rd
	            continue
	        elseif As_sw <= As_sw_rd
	            shear_reinforcements(row,:) = [stirrup_dia, t_spacing, As_sw_rd];
	            row=row+1;
	        end
	    end
	end
end

%% shear reinforcement
[~,order] = sort(shear_reinforcements(:,3));       
shear_reinforcements = shear_reinforcements(order,:);
shear_reinforcements_chosen = shear_reinforcements(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Torsional Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t, Ak, uk
tef = (widths*depths)/(2*(widths+depths));
Ak = (widths-tef)*(depths-tef);
uk = 2*(widths+depths-2*tef);

% interaction check
vv1 = 0.6*(1-fck/250);
T_rd_max = (1.33*vv1*fck*tef*Ak)/(cot_theta+(1/cot_theta));
check_T_V = (T_load/T_rd_max)+(V_load/V_rd_max);
if check_T_V > 1
	error('increase dimension of the beam to increase T_rd_max or V_rd_max')
end

% additional torsional stirrup reinforcement
As_st = (1e06)*T_load/(2*Ak*0.87*fyk*cot_theta);
st_max = min([(uk/8),0.75*d,min(widths,depths)]);

%% Overall Stirrup Design
As_s = As_sw+2*As_st;
s_max = min(sw_max,st_max);
% design spacing required
stirrup_reinforcements = [];
row = 1;
for stirrup_diameter = [8, 10, 12, 16, 20, 25, 32, 40]
    for t_spacing = s_max:-25:50
        As_s_rd = (2*pi*0.25*stirrup_diameter^2)/t_spacing;
        if As_s > As_s_rd
            continue
        elseif As_s <= As_s_rd
            stirrup_reinforcements(row,:) = [stirrup_diameter, t_spacing, As_s_rd];
            row=row+1;
        end
    end
end

%% stirrup reinforcement
[~,order] = sort(stirrup_reinforcements(:,3));       
stirrup_reinforcements = stirrup_reinforcements(order,:);
stirrup_reinforcements_chosen = stirrup_reinforcements(1,:);

%% additional Longitudinal Bars
As1_T = (1e06)*(T_load*uk*cot_theta)/(2*Ak*0.87*fyk);
rebar_areas_tor = [];
row_tor = 1;
for diameter_rebar = Rebar_dia
	As_per_bar = 0.25*pi*diameter_rebar^2;
	n_tor = ceil(As1_T/As_per_bar);
	A_tor_rd = As_per_bar*n_tor;
    if n_tor <= N_max_tor
        rebar_areas_tor(row_tor,:) = [As1_T, diameter_rebar, As_per_bar, n_tor, A_tor_rd];
        row_tor = row_tor+1;
    end
end

%% longitudinal reinforcement
[~,order] = sort(rebar_areas_tor(:,5));       
rebar_areas_tor = rebar_areas_tor(order,:);
rebar_areas_tor_chosen = rebar_areas_tor(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('(M+) longitudinal diameter = ')
disp(rebar_areas_pos_chosen)
disp('(M-) longitudinal diameter = ')
disp(rebar_areas_neg_chosen)
disp('(T) longitudinal diameter = ')
disp(rebar_areas_tor_chosen)
disp('(V + T) Stirrup Reinforcement = ')
disp(stirrup_reinforcements_chosen)