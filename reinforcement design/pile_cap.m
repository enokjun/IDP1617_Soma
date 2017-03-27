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
Date: 17th Feb 2017
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
%% soil
soil = dlmread('AA1_soil_profile.csv',',',0,0);

%% loads - unfactored
% for axial (P): compression (+) and tension (-)
% for moment (M): clockwise (+) and anti-clockwise (-)
col_number = 1;
loads = dlmread('hammersmith_core_loads.csv',',',0,0);
Dead_P = loads(1,col_number);         % unit: kN
Dead_Mx = loads(2,col_number);           % unit: kNm
Dead_My = loads(3,col_number);            % unit: kNm
Live_P = loads(4,col_number);         % unit: kN
Live_Mx = loads(5,col_number);           % unit: kNm
Live_My = loads(6,col_number);           % unit: kNm
Wind_P = loads(7,col_number);          % unit: kN
Wind_Mx = loads(8,col_number);           % unit: kNm
Wind_My = loads(9,col_number);            % unit: kNm

critical_P_load_combi = 49;
critical_Mx_load_combi = 49;
critical_My_load_combi = 29;
load_combination_case = [critical_P_load_combi, critical_Mx_load_combi, critical_My_load_combi];

%% factor of loads
Dead_A_Unf = 1.35;
Dead_A_F = 1;
Live_A_Unf = 1.5;
Wind_A_Unf = 0.75;

%% pile dimension
gamma_conc = 25;		% unit weight    unit: kN/m^3
d = 2.4;       % pile diameter     unit: m
N = 6;
Z = 45;		% pile depth    unit: m
L = 40;		% pile length    unit: m
dist_bw_piles = 2.5;
A = 0.25*pi*d^2;

%% col dimension
% col_width = 0.6;	 % unit: m
% col_length = 0.8;	 % unit: m

%% core dimension
core_area = 6.6;		% unit: m^2
ui = 53.6;				% punching shear perimeter   unit: m 

%% factors
gamma_c = 1.5;
gamma_s = 1.15;
alpha_cc = 0.85;

%% Concrete - C30/37 (EC2 TAle 3.1)
% strength 
fck = 30;        % unit: N/mm^2
fctm = 2.9;		 % unit: N/mm^2
cover = (25+75);      % unit: mm     % c_min + delta(c_dev)
rebar_space = 200;		% unit: mm
crack_control = 200;    % unit: mm  % crack control through limiting bar spaces
Ec = 33000;			 % unit: N/mm^2

%% Steel Rebar - S500 (EC2 TAle 3.1) 
% strength 
euk = 0.05;
fyk = 500;           % unit: N/mm^2
Es = 200000;        % unit: N/mm^2
% dimensions      
Rebar_dia = [12, 16, 20, 25, 32, 40];     % unit: mm     % bar diameter - minimum  
stirrup_diameter = [12, 16, 20, 25]; %[32,40];   % unit: mm

%% analysis
K0 = 0.168;

%% settlement limit
sett_limit = 1/1000;

%% cost ratio  
cost_rebar = 11775;		% unit: pound/m^3 -> based on 1500/tonne * 7.85 tonne/m^3
cost_rc = 200;		% unit: pound/m^3

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input from the ULS analysis - N 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assumption: centre-to-centre distance between piles (dx) is factor (2~3) x diameter
dx = dist_bw_piles*d;

switch N
%{
    case 4 % 2x2
        xn = dx*[-0.5; 0.5; -0.5; 0.5];
        yn = dx*[0.5; 0.5; -0.5; -0.5];
        pile_types = [1; 1; 1; 1];  
        which_cell_no = 1;
    case 5 % 2-1-2
        xn = dx*[-1; -1; 0; 1; 1];
        yn = dx*[1; -1; 0; 1; -1];
        pile_types = [1; 1; 2; 1; 1];  
        which_cell_no = 2;
%}
    case 6 % 2x3
        xn = dx*[-1; 0; 1; -1; 0; 1];
        yn = dx*[0.5; 0.5; 0.5; -0.5; -0.5; -0.5];
        pile_types = [1; 2; 1; 1; 2; 1];     % 1 = corner, 2 = edge
        % which_cell_no = 3;
%{
    case 7  %2-3-2
    	xn = dx*[-1; -1/2; -1/2; 0; 1/2; 1/2; 1];
        yn = dx*[0; -0.5*sqrt(3)*dx; 0.5*sqrt(3)*dx; 0; -0.5*sqrt(3)*dx; 0.5*sqrt(3)*dx; 0];x4
        pile_types = [1; 1; 1; 2; 1; 1; 1];     % 1 = corner, 2 = edge
        which_cell_no = 4;
    case 8 % 2x4
        xn = dx*[-1.5; -0.5; 0.5; 1.5;-1.5; -0.5; 0.5; 1.5];
        yn = dx*[0.5; 0.5; 0.5; 0.5; -0.5; -0.5; -0.5; -0.5];
        pile_types = [1; 2; 2; 1; 1; 2; 2; 1];     % 1 = corner, 2 = edge
        which_cell_no = 5;
    case 9  % 3x3 layout
        xn = dx*[-1; 0; 1; -1; 0; 1; -1; 0; 1];
        yn = dx*[1; 1; 1; 0; 0; 0; -1; -1; -1];
        pile_types = [1; 2; 1; 2; 3; 2; 1; 2; 1];
        which_cell_no = 6; 
%}
    otherwise
        error('something wrong too many piles or not enough - check your codes')
end

%% Ixx and Iyy
Ixx = 0;
Iyy = 0;
for qqqq = 1:length(xn)
    Ixx = Ixx + yn(qqqq)^2;
    Iyy = Iyy + xn(qqqq)^2;
end

%% group
numberofpilesinx = length(unique(xn));
numberofpilesiny = length(unique(yn));
bb = 0.15*2+d+(numberofpilesinx-1)*dx;
ll = 0.15*2+d+(numberofpilesiny-1)*dx;     

%% distribution of axial loads
pile_coordinates = [xn, yn];       % coordinate showing pile locations
N_pile = [];
unique_piles = unique(pile_types);
for unique_piles_analyse = 1:length(unique_piles)
    N_pile(unique_piles_analyse) = [sum(pile_types == unique_piles(unique_piles_analyse))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load combination analysis - total loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% factored loads
ULS_Dead_Unf = Dead_A_Unf*[Dead_P, Dead_Mx, Dead_My];
ULS_Dead_F = Dead_A_F*[Dead_P, Dead_Mx, Dead_My];
ULS_Live_Unf = Live_A_Unf*[Live_P, Live_Mx, Live_My];
ULS_Wind_Unf = Wind_A_Unf*[Wind_P, Wind_Mx, Wind_My];
SLS_Dead_Unf = [Dead_P, Dead_Mx, Dead_My];
SLS_Live_Unf = [Live_P, Live_Mx, Live_My];
SLS_Wind_Unf = [Wind_P, Wind_Mx, Wind_My];

%% find the critical load combinations
ULS_critical_loads = [];
SLS_critical_loads = [];
for combi_PMM = 1:3

	critical_load_combinations = load_combination_case(combi_PMM);

	switch critical_load_combinations
	    case 19
	    	ULS_Load_comp_total = ULS_Dead_Unf(1);
	        % ULS_Load_ten_total = 0;
	        ULS_Load_Mx = ULS_Dead_Unf(2);
	        ULS_Load_My = ULS_Dead_Unf(3);
	        SLS_Load_comp_total = SLS_Dead_Unf(1);
	        % SLS_Load_ten_total = 0;
	        SLS_Load_Mx = SLS_Dead_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3);
	    case 91
	    	ULS_Load_comp_total = 0;
	        % ULS_Load_ten_total = ULS_Dead_Unf(1); 
	        ULS_Load_Mx = ULS_Dead_Unf(2);
	        ULS_Load_My = ULS_Dead_Unf(3);
	        SLS_Load_comp_total = 0;
	        % SLS_Load_ten_total = SLS_Dead_Unf(1); 
	        SLS_Load_Mx = SLS_Dead_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3);
	    case 13
	    	ULS_Load_comp_total = ULS_Dead_Unf(1); 
	        % ULS_Load_ten_total = ULS_Dead_F(1)+ULS_Live_Unf(1); 
	        ULS_Load_Mx = max(ULS_Dead_Unf(2), ULS_Dead_F(2)+ULS_Live_Unf(2));
	        ULS_Load_My = max(ULS_Dead_Unf(3), ULS_Dead_F(3)+ULS_Live_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_Unf(1); 
	        % SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	    case 31
	    	ULS_Load_comp_total = ULS_Dead_F(1)+ULS_Live_Unf(1); 
	        % ULS_Load_ten_total = ULS_Dead_Unf(1); 
	        ULS_Load_Mx = max(ULS_Dead_Unf(2), ULS_Dead_F(2)+ULS_Live_Unf(2));
	        ULS_Load_My = max(ULS_Dead_Unf(3), ULS_Dead_F(3)+ULS_Live_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        % SLS_Load_ten_total = SLS_Dead_Unf(1); 
	        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	    case 15
	    	ULS_Load_comp_total = ULS_Dead_Unf(1); 
	        % ULS_Load_ten_total = ULS_Dead_F(1)+ULS_Live_Unf(1)+ULS_Wind_Unf(1); 
	        ULS_Load_Mx = max(ULS_Dead_Unf(2), ULS_Dead_F(2)+ULS_Live_Unf(2)+ULS_Wind_Unf(2));
	        ULS_Load_My = max(ULS_Dead_Unf(3), ULS_Dead_F(3)+ULS_Live_Unf(3)+ULS_Wind_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_Unf(1); 
	        % SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
	        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3));
	    case 51
	    	ULS_Load_comp_total = ULS_Dead_F(1)+ULS_Live_Unf(1)+ULS_Wind_Unf(1); 
	        % ULS_Load_ten_total = ULS_Dead_Unf(1); 
	        ULS_Load_Mx = max(ULS_Dead_Unf(2), ULS_Dead_F(2)+ULS_Live_Unf(2)+ULS_Wind_Unf(2));
	        ULS_Load_My = max(ULS_Dead_Unf(3), ULS_Dead_F(3)+ULS_Live_Unf(3)+ULS_Wind_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
	        % SLS_Load_ten_total = SLS_Dead_Unf(1); 
	        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3));
	    case 29
	    	ULS_Load_comp_total = ULS_Dead_Unf(1)+ULS_Live_Unf(1);
	        % ULS_Load_ten_total = 0;
	        ULS_Load_Mx = ULS_Dead_Unf(2)+ULS_Live_Unf(2);
	        ULS_Load_My = ULS_Dead_Unf(3)+ULS_Live_Unf(3);
	        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
	        % SLS_Load_ten_total = 0;
	        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3);
	    case 92
	    	ULS_Load_comp_total = 0;
	        % ULS_Load_ten_total = ULS_Dead_Unf(1)+ULS_Live_Unf(1);
	        ULS_Load_Mx = ULS_Dead_Unf(2)+ULS_Live_Unf(2);
	        ULS_Load_My = ULS_Dead_Unf(3)+ULS_Live_Unf(3);
	        SLS_Load_comp_total = 0;
	        % SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
	        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3);
	    case 49
	    	ULS_Load_comp_total = ULS_Dead_Unf(1)+ULS_Live_Unf(1)+ULS_Wind_Unf(1);
	        % ULS_Load_ten_total = 0;
	        ULS_Load_Mx = ULS_Dead_Unf(2)+ULS_Live_Unf(2)+ULS_Wind_Unf(2);
	        ULS_Load_My = ULS_Dead_Unf(3)+ULS_Live_Unf(3)+ULS_Wind_Unf(3);
	        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
	        % SLS_Load_ten_total = 0;
	        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3);
	    case 94
	    	ULS_Load_comp_total = 0;
	        % ULS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
	        ULS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2);
	        ULS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3);
	        SLS_Load_comp_total = 0;
	        % SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
	        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2);
	        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3);
	    case 36
	    	ULS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        % ULS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
	        ULS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        ULS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        % SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
	        SLS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	    case 63
	    	ULS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
	        % ULS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        ULS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        ULS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
	        % SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
	        SLS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
	        SLS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
	end

	if combi_PMM == 1
		ULS_critical_loads(1) = ULS_Load_comp_total;
		SLS_critical_loads(1) = SLS_Load_comp_total;
	elseif combi_PMM == 2
		ULS_critical_loads(2) = ULS_Load_Mx;
		SLS_critical_loads(2) = SLS_Load_Mx;
	elseif combi_PMM == 3
		ULS_critical_loads(3) = ULS_Load_My;
		SLS_critical_loads(3) = SLS_Load_My;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ULS additional load on each piles based on moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Mx
if ULS_critical_loads(2) > 0
	Mx_factor = 1;
elseif ULS_critical_loads(2) < 0
	Mx_factor = -1;
elseif ULS_critical_loads(2) == 0
	Mx_factor = 0;
end
% My
if ULS_critical_loads(3) > 0
	My_factor = 1;
elseif ULS_critical_loads(3) < 0
	My_factor = -1;
elseif ULS_critical_loads(3) == 0
	My_factor = 0;
end

ULS_loads_from_M = []; 
for m = 1:length(pile_types)	

	if xn(m) > 0
		xn_factor = 1;
	elseif xn(m) < 0
		xn_factor = -1;
	elseif xn(m) == 0
		xn_factor = 0;
	end
	if xn(m) > 0
		yn_factor = 1;
	elseif xn(m) < 0
		yn_factor = -1;
	elseif xn(m) == 0
		yn_factor = 0;
	end

    ULS_loads_from_M(m) = Mx_factor*yn_factor*(ULS_critical_loads(2)*yn(m))/Ixx + My_factor*xn_factor*(ULS_critical_loads(3)*xn(m))/Iyy;                      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLS additional load on each piles based on moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Mx
if SLS_critical_loads(2) > 0
	Mx_factor = 1;
elseif SLS_critical_loads(2) < 0
	Mx_factor = -1;
elseif SLS_critical_loads(2) == 0
	Mx_factor = 0;
end
% My
if SLS_critical_loads(3) > 0
	My_factor = 1;
elseif SLS_critical_loads(3) < 0
	My_factor = -1;
elseif SLS_critical_loads(3) == 0
	My_factor = 0;
end

SLS_loads_from_M = []; 
for m = 1:length(pile_types)	

	if xn(m) > 0
		xn_factor = 1;
	elseif xn(m) < 0
		xn_factor = -1;
	elseif xn(m) == 0
		xn_factor = 0;
	end
	if xn(m) > 0
		yn_factor = 1;
	elseif xn(m) < 0
		yn_factor = -1;
	elseif xn(m) == 0
		yn_factor = 0;
	end

    SLS_loads_from_M(m) = Mx_factor*yn_factor*(SLS_critical_loads(2)*yn(m))/Ixx + My_factor*xn_factor*(SLS_critical_loads(3)*xn(m))/Iyy;                      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ULS - load distribution from group pile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% soil Young's and Shear Modulus
if Z >= soil(6,1) && Z < soil(6,1)+soil(6,2) 
    Eb = soil(6,8)*1000;
    v = soil(6,9);
elseif Z >= soil(5,1) && Z < soil(6,1) 
    Eb = soil(5,8)*1000;
    v = soil(5,9);
elseif Z >= soil(4,1) && Z < soil(5,1) 
    Eb = soil(4,8)*1000;
    v = soil(4,9);
elseif Z >= soil(3,1) && Z < soil(4,1) 
    Eb = soil(3,8)*1000;
    v = soil(3,9);
elseif Z >= soil(2,1) && Z < soil(3,1) 
    Eb = soil(2,8)*1000;
    v = soil(2,9);
elseif Z >= soil(1,1) && Z < soil(2,1) 
    Eb = soil(1,8)*1000;   
    v = soil(1,9);
end
G = Eb/(2*(1+v));

%% Fa
mu = sqrt((2*pi*G)/(log(2*L/d)*Eb*A));
Kbi = (2*d*G)/(1-v);
omega = Kbi/(Eb*A*mu);
mu_L_2 = 2*mu*L; 
Fa = (mu_L_2 + sinh(mu_L_2)+(omega^2)*(sinh(mu_L_2)-mu_L_2)+2*omega*(cosh(mu_L_2)-1))/((2+2*omega^2)*sinh(mu_L_2)+4*omega*cosh(mu_L_2));

%% aj
list_of_piles_to_calc = unique(pile_types);

aj_per_eq = [];
k = 1;

for i_loops = 1:length(list_of_piles_to_calc)  % shows how many pile types are requried 
    % calculate aj for each pile types
    aj_sum_each_pile_type = zeros(1,length(N_pile));
    for m = 1:length(pile_types)
        two_points = [pile_coordinates(list_of_piles_to_calc(i_loops),:); pile_coordinates(m,:)];
        s = pdist(two_points);
        if s == 0
            aj = 1;
        else
            aj = Fa*(log(L/s)/log(2*L/d));
        end
        aj_sum_each_pile_type(pile_types(m)) = aj_sum_each_pile_type(pile_types(m)) + aj;
    end
    aj_per_eq(k,:) = aj_sum_each_pile_type;
    k=k+1;
end

% assume rigid pile cap; therefore, all settlement are equal
number_of_pile_types_to_analyse = length(nonzeros(unique(pile_types)));

M = aj_per_eq;
% check matrix singularity - for now skip this in a loop
if det(M) == 0
    error('something wrong with the matrix M')
end

%ratio_of_loads = inv(M)*ones(number_of_pile_types_to_analyse,1);
ratio_of_loads = M\ones(number_of_pile_types_to_analyse,1);
Q1_ULS = ULS_critical_loads(1)/(N_pile*ratio_of_loads);
ULS_load_P_per_pile_type_P = Q1_ULS*ratio_of_loads;

Q1_SLS = SLS_critical_loads(1)/(N_pile*ratio_of_loads);
SLS_load_P_per_pile_type_P = Q1_SLS*ratio_of_loads;

% calculations check
sum_of_individual_loads = sum(N_pile*ULS_load_P_per_pile_type_P);
if round(ULS_critical_loads(1)) ~= round(sum_of_individual_loads)
    error('your sum of loads on each pile does not add up to total uls axial loads you are imposing')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loads on each pile types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ULS_load_on_each_pile = [];
SLS_load_on_each_pile = [];
for piles_to_loop = 1:length(ULS_loads_from_M) 
	ULS_load_on_each_pile(piles_to_loop) = ULS_loads_from_M(piles_to_loop) + ULS_load_P_per_pile_type_P(pile_types(piles_to_loop));
	SLS_load_on_each_pile(piles_to_loop) = SLS_loads_from_M(piles_to_loop) + SLS_load_P_per_pile_type_P(pile_types(piles_to_loop));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design to find the pile cap thickness required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design_based_t = [];
new_line = 1;

for t = 1:0.1:3

	%% pile cap 
	pile_cap_V = bb*ll*t;
	pile_cap_W = gamma_conc*pile_cap_V;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Moment design X-X - write for pile = 6, generalise later
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Mxx_pile = (ll/2)*max([ULS_load_on_each_pile(1)+ULS_load_on_each_pile(4); ULS_load_on_each_pile(3)+ULS_load_on_each_pile(6)]);
	Mxx_pile_cap = 0.125*(pile_cap_W/ll)*(ll-d)^2;
	Mxx = Mxx_pile - Mxx_pile_cap;

	% assume two layers of rebar diameter 12mm at the start
	d_pos = 1000*t - (Rebar_dia(1) + stirrup_diameter(1) + cover);
	d_prime_pos = Rebar_dia(1) + stirrup_diameter(1) + cover;
	K_pos = (1e06)*Mxx/(fck*(bb*1000)*d_pos^2);
	if K_pos <= K0
		As_pos_xx = (1e06)*Mxx/(0.87*fyk*0.8*d_pos);
		As_neg_xx = 0;
	elseif K_pos > K0
		M_rd_max = (1e-06)*K0*fck*(bb*1000)*d_pos^2;
		if d_prime_pos/(0.45*d_pos) > 0.38
			continue
			% error('decrease d_prime_pos')
		end
		As_neg_xx = (1e06)*(Mxx-M_rd_max)/(0.87*fyk*(d_pos-d_prime_pos));
		ZZ = d_pos*(0.5+sqrt(0.25-(3/3.4)*K0));
		As_pos_xx = (1e06)*(Mxx/(0.87*fyk*ZZ))+As_neg_xx;
	end

	%% Required Postive and negative reinforcement areas
	rebars_to_design_XX = [As_pos_xx, As_neg_xx];
	rebar_areas_pos_XX = [];
	rebar_areas_neg_XX  = [];
	row_p = 1;
	row_n = 1;
	
	for designing_long_As = 1:2
		for diameter_rebar_1 = Rebar_dia
			for diameter_rebar_2 = Rebar_dia
				As_per_bar_1 = 0.25*pi*diameter_rebar_1^2;
				As_per_bar_2 = 0.25*pi*diameter_rebar_2^2;

				l_spacing_max = min(250,crack_control);
				l_spacing_min_1 = max([diameter_rebar_1,25,20]);		% assume: min aggregate size = 20mm
				l_spacing_min_2 = max([diameter_rebar_2,25,20]);

				for s_line1 = l_spacing_max:-50:l_spacing_min_1
					for s_line2 = l_spacing_max:-50:l_spacing_min_2
						N_line1 = floor(bb*1000/s_line1);
						N_line2 = floor(bb*1000/s_line2);

						As_rd_xx = As_per_bar_1*N_line1 + As_per_bar_2*N_line2;

						if rebars_to_design_XX(designing_long_As) <= As_rd_xx
							if designing_long_As == 1
								rebar_areas_pos_XX(row_p,:) = [As_pos_xx, diameter_rebar_1, s_line1, N_line1, diameter_rebar_2, s_line2, N_line2, As_rd_xx];
						        row_p = row_p+1;
						    elseif designing_long_As == 2
								rebar_areas_neg_XX(row_n,:) = [As_neg_xx, diameter_rebar_1, s_line1, N_line1, diameter_rebar_2, s_line2, N_line2, As_rd_xx];
						        row_n = row_n+1;
					       	end
					    else 
					    	continue
					    end
					end
				end
		    end
		end
	end

    if isempty(rebar_areas_pos_XX) == 1 ||  isempty(rebar_areas_neg_XX) == 1
        continue
    end
    
	[~,order] = sort(rebar_areas_pos_XX(:,8));       
	rebar_areas_pos_XX  = rebar_areas_pos_XX(order,:);
	[~,order] = sort(rebar_areas_neg_XX(:,8));       
	rebar_areas_neg_XX  = rebar_areas_neg_XX(order,:);

	% final reinforcement area
	rebar_areas_pos_chosen_XX = rebar_areas_pos_XX(1,:);
	rebar_areas_neg_chosen_XX = rebar_areas_neg_XX(1,:);

	%% Check - longitudinal area
	As_min = max(0.26*fctm*(bb*1000)*d_pos/fyk, 0.0015*(bb*1000)*d_pos);
	As_max = 0.04*(bb*1000)*(t*1000);

	if rebar_areas_pos_chosen_XX(8) < As_min 
		new_As = find(rebar_areas_pos_XX(:,8) >= As_min);
		rebar_areas_pos_chosen_XX = rebar_areas_pos_XX(new_As(1),:);
	elseif rebar_areas_pos_chosen_XX(8) > As_max
		new_As = find(rebar_areas_pos_XX(:,8) <= As_max);
		rebar_areas_pos_chosen_XX = rebar_areas_pos_XX(new_As(1),:);
	end

	if rebar_areas_neg_chosen_XX(8) < As_min 
		new_As = find(rebar_areas_neg_XX(:,8) >= As_min);
		rebar_areas_neg_chosen_XX = rebar_areas_neg_XX(new_As(1),:);
	elseif rebar_areas_neg_chosen_XX(8) > As_max
		new_As = find(rebar_areas_neg_XX(:,8) <= As_max);
		rebar_areas_neg_chosen_XX = rebar_areas_neg_XX(new_As(1),:);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Shear design X-X - write for pile = 6, generalise later
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% critical face is 20% inward from the pile diameter
	Vxx_pile = max([ULS_load_on_each_pile(1)+ULS_load_on_each_pile(4); ULS_load_on_each_pile(3)+ULS_load_on_each_pile(6)]);
	Vxx_pile_cap = (pile_cap_W/ll)*(0.8*d + 0.15);
	Vxx = Vxx_pile - Vxx_pile_cap; 

	%% Shear Resistance Reinforcement
	ddx = (t*1000)-cover-stirrup_diameter(1)-rebar_areas_pos_chosen_XX(2);
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
	z = 0.9*ddx;
	% V_rd_max
	V_rd_max = (bb*v1*fcd*z)/(cot_theta+(1/cot_theta));
	if Vxx > V_rd_max
	    continue
	    % error('try a different diameter or something... Vxx > V_rd_max')
	end

	% As/sw
	As_sw_xx = 1000*Vxx/(z*fyd*cot_theta);
	% sw_max - round down to closest 25mm
	sw_max = floor((0.75*ddx)/25)*25;

	%% design spacing required
	shear_reinforcements_xx = [];
	row = 1;
	for stirrup_dia = stirrup_diameter
	    for t_spacing = sw_max:-25:50
	        As_sw_rd = (2*pi*0.25*stirrup_dia^2)/t_spacing;
	        if As_sw_xx > As_sw_rd
	            continue
	        elseif As_sw_xx <= As_sw_rd
	            shear_reinforcements_xx(row,:) = [stirrup_dia, t_spacing, As_sw_rd];
	            row=row+1;
	        end
	    end
    end
    
    if isempty(shear_reinforcements_xx) == 1
        continue
    end
    
	%% shear reinforcement
	[~,order] = sort(shear_reinforcements_xx(:,3));       
	shear_reinforcements_xx = shear_reinforcements_xx(order,:);
	shear_reinforcements_chosen_xx = shear_reinforcements_xx(1,:);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Moment design Y-Y - write for pile = 6, generalise later
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Myy_pile = (bb/2)*max([ULS_load_on_each_pile(1)+ULS_load_on_each_pile(2)+ULS_load_on_each_pile(3); ULS_load_on_each_pile(4)+ULS_load_on_each_pile(5)+ULS_load_on_each_pile(6)]);
	Myy_pile_cap = 0.125*(pile_cap_W/bb)*(bb-d)^2;
	Myy = Myy_pile - Myy_pile_cap;

	% assume two layers of rebar diameter 12mm at the start
	d_pos = 1000*t - (Rebar_dia(1) + stirrup_diameter(1) + cover);
	d_prime_pos = Rebar_dia(1) + stirrup_diameter(1) + cover;
	K_pos = (1e06)*Myy/(fck*(ll*1000)*d_pos^2);
	if K_pos <= K0
		As_pos_yy = (1e06)*Myy/(0.87*fyk*0.8*d_pos);
		As_neg_yy = 0;
	elseif K_pos > K0
		M_rd_max = (1e-06)*K0*fck*(ll*1000)*d_pos^2;
		if d_prime_pos/(0.45*d_pos) > 0.38
			continue
			% error('decrease d_prime_pos')
		end
		As_neg_yy = (1e06)*(Myy-M_rd_max)/(0.87*fyk*(d_pos-d_prime_pos));
		ZZ= d_pos*(0.5+sqrt(0.25-(3/3.4)*K0));
		As_pos_yy = (1e06)*(Myy/(0.87*fyk*ZZ))+As_neg_yy;
	end

	%% Required Postive and negative reinforcement areas
	rebars_to_design_YY = [As_pos_yy, As_neg_yy];
	rebar_areas_pos_YY = [];
	rebar_areas_neg_YY  = [];
	row_p = 1;
	row_n = 1;

    for designing_long_As = 1:2
		for diameter_rebar_1 = Rebar_dia
			for diameter_rebar_2 = Rebar_dia
				As_per_bar_1 = 0.25*pi*diameter_rebar_1^2;
				As_per_bar_2 = 0.25*pi*diameter_rebar_2^2;

				l_spacing_max = min(250,crack_control);
				l_spacing_min_1 = min([diameter_rebar_1,25,20]);		% assume: min aggregate size = 20mm
				l_spacing_min_2 = min([diameter_rebar_2,25,20]);	

				for s_line1 = l_spacing_max:-50:l_spacing_min_1
					for s_line2 = l_spacing_max:-50:l_spacing_min_2
						N_line1 = floor(ll*1000/s_line1);
						N_line2 = floor(ll*1000/s_line2);
						As_rd_yy = As_per_bar_1*N_line1 + As_per_bar_2*N_line2;

						if rebars_to_design_YY(designing_long_As) <= As_rd_yy
							if designing_long_As == 1
								rebar_areas_pos_YY(row_p,:) = [As_pos_yy, diameter_rebar_1, s_line1, N_line1, diameter_rebar_2, s_line2, N_line2, As_rd_yy];
						        row_p = row_p+1;
						    elseif designing_long_As == 2
								rebar_areas_neg_YY(row_n,:) = [As_neg_yy, diameter_rebar_1, s_line1, N_line1, diameter_rebar_2, s_line2, N_line2, As_rd_yy];
						        row_n = row_n+1;
					       	end
					    else 
					    	continue
					    end
					end
				end
		    end
		end
	end

    if isempty(rebar_areas_pos_YY) == 1 ||  isempty(rebar_areas_neg_YY) == 1
        continue
    end
    
	[~,order] = sort(rebar_areas_pos_YY(:,8));       
	rebar_areas_pos_YY = rebar_areas_pos_YY(order,:);
	[~,order] = sort(rebar_areas_neg_YY(:,8));       
	rebar_areas_neg_YY = rebar_areas_neg_YY(order,:);

	% final reinforcement area
	rebar_areas_pos_chosen_YY = rebar_areas_pos_YY(1,:);
	rebar_areas_neg_chosen_YY = rebar_areas_neg_YY(1,:);

	%% Check - longitudinal area
	As_min = max(0.26*fctm*(ll*1000)*d_pos/fyk, 0.0015*(ll*1000)*d_pos);
	As_max = 0.04*(ll*1000)*(t*1000);

	if rebar_areas_pos_chosen_YY(8) < As_min 
		new_As = find(rebar_areas_pos_YY(:,8) >= As_min);
		rebar_areas_pos_chosen_YY = rebar_areas_pos_YY(new_As(1),:);
	elseif rebar_areas_pos_chosen_YY(8) > As_max
		new_As = find(rebar_areas_pos_YY(:,8) <= As_max);
		rebar_areas_pos_chosen_YY = rebar_areas_pos_YY(new_As(1),:);
	end

	if rebar_areas_neg_chosen_YY(8) < As_min 
		new_As = find(rebar_areas_neg_YY(:,8) >= As_min);
		rebar_areas_neg_chosen_YY = rebar_areas_neg_YY(new_As(1),:);
	elseif rebar_areas_neg_chosen_YY(8) > As_max
		new_As = find(rebar_areas_neg_YY(:,8) <= As_max);
		rebar_areas_neg_chosen_YY = rebar_areas_neg_YY(new_As(1),:);
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Shear design Y-Y - write for pile = 6, generalise later
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% critical face is 20% inward from the pile diameter
	Vyy_pile = max([ULS_load_on_each_pile(1)+ULS_load_on_each_pile(2)+ULS_load_on_each_pile(3); ULS_load_on_each_pile(4)+ULS_load_on_each_pile(5)+ULS_load_on_each_pile(6)]);
	Vyy_pile_cap = (pile_cap_W/bb)*(0.8*d + 0.15);
	Vyy = Vyy_pile - Vyy_pile_cap; 

	%% Shear Resistance Reinforcement
	ddy = (t*1000)-cover-stirrup_diameter(1)-rebar_areas_pos_chosen_YY(2);
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
	z = 0.9*ddy;
	% V_rd_max
	V_rd_max = (ll*v1*fcd*z)/(cot_theta+(1/cot_theta));
	if Vyy > V_rd_max
	    % disp('try a different diameter or something... Vyy > V_rd_max')
	    continue
	end

	% As/sw
	As_sw_yy = 1000*Vyy/(z*fyd*cot_theta);
	% sw_max - round down to closest 25mm
	sw_max = floor(0.75*ddy/25)*25;

	%% design spacing required
	shear_reinforcements_yy = [];
	row = 1;
	for stirrup_dia = stirrup_diameter
	    for t_spacing = sw_max:-25:50
	        As_sw_rd = (2*pi*0.25*stirrup_dia^2)/t_spacing;
	        if As_sw_yy > As_sw_rd
	            continue
	        elseif As_sw_yy <= As_sw_rd
	            shear_reinforcements_yy(row,:) = [stirrup_dia, t_spacing, As_sw_rd];
	            row=row+1;
	        end
	    end
    end
    
    if isempty(shear_reinforcements_yy) == 1
        continue
    end

	%% shear reinforcement
	[~,order] = sort(shear_reinforcements_yy(:,3));       
	shear_reinforcements_yy = shear_reinforcements_yy(order,:);
	shear_reinforcements_chosen_yy = shear_reinforcements_yy(1,:);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Punching Shear design - write for pile = 6, generalise later
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dd = (ddx+ddy)/2;

	% ved_col = 
	ved_core = 1.15*ULS_critical_loads(1)/(ui*dd);
	ved_pile = 1.5*max(ULS_load_on_each_pile)/(ui*dd);
	ved = max(ved_core,ved_pile);

	%% concrete shear design
	Crd_c = 0.18/gamma_c;	
	k = min(1+sqrt(200/dd),2);		
	rho1 = min(sqrt((rebar_areas_neg_chosen_XX(5)/(bb*dd))+(rebar_areas_neg_chosen_YY(5)/(ll*dd))),0.02);	
	vmin = 0.035*(k^1.5)*sqrt(fck);
	sigma_cp = min(ULS_critical_loads(1)/core_area, 0.2*fcd);
	vrd_c = max(Crd_c*k*((100*rho1*fck)^(1/3))+0.1*sigma_cp, vmin+0.1*sigma_cp);

	% vrd_max = min(0.9*sqrt(fck),0.9);
	vrd_max = min(2*vrd_c,0.5*(0.6*(1-(fck/250)))*fcd);

	fyd_eff = min(250+0.25*dd,fyk);

	% Does it require additional punching steel reinforcement?
	if vrd_c >= ved
		% disp('no need for additional shear reinforcement from the punching shear')
        punching_shear_reinforcements_chosen = [0, 0, 0];
	elseif vrd_c < ved && vrd_max >= ved
		% As_sw
		As_sr = (ved-0.75*vrd_c)/(u1*1.5*fyd_eff);
		% sw_max - round down to closest 25mm
		sr_max = floor(0.75*dd/25)*25;

		%% design spacing required
		punching_shear_reinforcements = [];
		row = 1;
		for stirrup_dia = stirrup_diameter
		    for t_spacing = sr_max:-25:50
		        As_sr_rd = (2*pi*0.25*stirrup_dia^2)/t_spacing;
		        if As_sr > As_sw_rd
		            continue
		        elseif As_sr <= As_sr_rd
		            punching_shear_reinforcements(row,:) = [stirrup_dia, t_spacing, As_sr_rd];
		            row=row+1;
		        end
		    end
		end
		%% shear reinforcement
		[~,order] = sort(punching_shear_reinforcements(:,3));       
		punching_shear_reinforcements = punching_shear_reinforcements(order,:);
		punching_shear_reinforcements_chosen = punching_shear_reinforcements(1,:);
	elseif vrd_max < ved
		continue
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% deflection limit
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Ixx_t = bb*(t^3)/12;			% unit: m^4
	Iyy_t = ll*(t^3)/12;			% unit: m^4
	delta_xx = 0.001*SLS_critical_loads(1)*(dx^3)/(48*Ec*1000*Ixx);
	delta_yy = 0.001*SLS_critical_loads(1)*((2*dx)^3)/(48*Ec*1000*Iyy);
	delta_sum = delta_yy+delta_xx;

	if delta_sum/(0.5*ll*1000) > sett_limit || delta_sum/(0.5*bb*1000) > sett_limit 
		continue
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% calculate cost of pile cap depending - concrete and steel volume
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% rebar cross-sectional area used  unit: mm^3
	long_V_total = 1000*bb*(rebar_areas_pos_chosen_XX(end)+rebar_areas_neg_chosen_XX(end))+1000*ll*(rebar_areas_pos_chosen_YY(end)++rebar_areas_neg_chosen_YY(end));
	stirrup_V_xx = 1000*bb*(shear_reinforcements_chosen_xx(3)*shear_reinforcements_chosen_xx(2)*(ll*1000/shear_reinforcements_chosen_xx(2)));
	stirrup_V_yy = 1000*ll*(shear_reinforcements_chosen_yy(3)*shear_reinforcements_chosen_yy(2)*(bb*1000/shear_reinforcements_chosen_yy(2)));
    if isequal(punching_shear_reinforcements_chosen,zeros(1,3))
        stirrup_V_punching_shear = 0;
    else
        stirrup_V_punching_shear = 1000*bb*(punching_shear_reinforcements_chosen(3)*punching_shear_reinforcements_chosen(2)*(ll*1000/punching_shear_reinforcements_chosen(2)))+...
							   1000*ll*(punching_shear_reinforcements_chosen(3)*punching_shear_reinforcements_chosen(2)*(bb*1000/punching_shear_reinforcements_chosen(2))); 
    end
    
	% rebar cost
	Long_V_cost = (1e-09)*long_V_total*cost_rebar;
	Stirrup_V_cost = (1e-09)*(stirrup_V_xx+stirrup_V_yy+stirrup_V_punching_shear)*cost_rebar;
	Overall_Steel_V_cost = Long_V_cost+Stirrup_V_cost;

	% RC pile cap cost
	RC_pile_cap_cost = pile_cap_V*cost_rc;

	% overall cost
	Overall_pile_cap_cost = RC_pile_cap_cost+Overall_Steel_V_cost;

	costs = [RC_pile_cap_cost, Long_V_cost, Stirrup_V_cost, Overall_Steel_V_cost, Overall_pile_cap_cost];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% tabulate design
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	design_based_t(new_line,:) = [bb, ll, t, rebar_areas_pos_chosen_XX, rebar_areas_pos_chosen_YY, rebar_areas_neg_chosen_XX, rebar_areas_neg_chosen_YY, shear_reinforcements_chosen_xx, shear_reinforcements_chosen_yy, punching_shear_reinforcements_chosen, costs];
	new_line = new_line + 1;
end

[~,order] = sort(design_based_t(:,end));       
design_based_t = design_based_t(order,:);

dlmwrite('pile_cap_design_Hammersmith_core.csv',design_based_t,'delimiter',',')
disp('done')