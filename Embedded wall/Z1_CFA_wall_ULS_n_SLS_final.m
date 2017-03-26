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
Date: 5th Dec 2016
Use: IDP CFA bored embedded wall capacity design
Key Assumptions:
	- concrete class C25/30
	- Using Design Approach 1 (EC7 2.4.7.3.4.2)
	- Fleming's Hypebolic Method for pile settlement
	- Thickness and depth of soil profile are generalised to be 
	  leveled, based on the borehole logs
    - Beta method for skin friction on fill and terrace gravel
        - Beta assumed as Ks =0.7 and 
                fill interface-friction-angle = 17
                sand interface-friction-angle = 22

References:
    - BS EN 1997:1-2004 (EC7)
    - BS NA EN 1997:1-2004 (BS ANNEX EC7)
    - IStructE Manual for the geotechnical design of structures to Eurocode 7
    - LDSA Foundation No 1: guidance notes for the design of straight 
        shafted bored piles in london clay - 2009 and 2000 ed.
    - Geotechnical Earthquake Engineering by S.L. Kramer 
    - CIRIA C580 Guide
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% column location names
column_list = {'XX 2', 'XX 3'};

% Design load
Dead_P_list = [-1000, 1000];      % unit: kN
Live_P_list = [-1000, 300];       % unit: kN
Wind_P_list = [-500, 1000];       % unit: kN

load_type_list = {'compression','tension'};

% soil
% [depth/m, thickness/m, bulk density/kNm^-3, Cu/kPa, friction angle, c/kPa, alpha, E/MPa, v, Cr, Cc]
soil_list = {'AA1_soil_profile.csv','AA2_soil_profile.csv'};

% basement dimension
basement_depth = 4;	% unit: m
soil_otherside_depth = 0;   % unit: m
column_length_span = 8; % unit: m
column_width_span = 10; % unit: m

% inital pile design range
B = [0.6, 0.75, 0.9, 1, 1.2];  % unit: m
Lmin = 10;   % unit: m
Lmax = 30;   % unit: m
Lint = 0.5;  % increment size

% CIRIA C580 Guide
Mmin = 3;
Mmax = 7; 
Space = [0.9, 1.15, 1.35, 1.5, 1.8];  % unit: m

%% STR/GEO Load combination factor (EC7 Table A.3)
% A1
A1_STR_GEO_Dead_Unf = 1.35;     % dead - unfavourable
A1_STR_GEO_Dead_F = 1.0;        % dead - favourable
A1_STR_GEO_Live = 1.5;          % if favourable = 0
A1_STR_GEO_wind = 1.5*0.5;      % if favourable = 0
% A2
A2_STR_GEO_Dead_Unf = 1.0;     % unfavourable
A2_STR_GEO_Dead_F = 1.0;       % favourable
A2_STR_GEO_Live = 1.3;         % if favourable = 0
A2_STR_GEO_wind = 1.3*0.5;     % if favourable = 0

% Material factors (EC7 TAle A.4)
% M1
M1_STR_GEO_phi = 1.0;   % friction angle
M1_STR_GEO_c = 1.0;     % cohesion
M1_STR_GEO_Cu = 1.0;    % undrained shear strength
M1_STR_GEO_den = 1.0;   % weight density
% M2
M2_STR_GEO_phi = 1.25;  % friction angle
M2_STR_GEO_c = 1.25;    % cohesion
M2_STR_GEO_Cu = 1.4;    % undrained shear strength
M2_STR_GEO_den = 1.4;   % weight density

% Partial Resistance factors (BS Annex EC7 Table A.NA.7)
% R1
R1_STR_GEO_base = 1.0;
R1_STR_GEO_sh_comp = 1.0;
R1_STR_GEO_sh_ten = 1.0;
R1_STR_GEO_comb_comp = 1.0;
% R4 - no explicit verification 
R4_STR_GEO_base = 2.0;
R4_STR_GEO_sh_comp = 1.6;
R4_STR_GEO_sh_ten = 2.0;
R4_STR_GEO_comb_comp = 2.0;

% beta method parameters 
Ks = 0.7;           %(IStructE Manuel for EC7, TAle 7.24) 
del_fill = 17;      %(Kramer, Table 11-1)
del_tergravel = 22; %(Kramer, Table 11-1)

% Concrete - C30/37 (EC2 TAle 3.1) 
fck_cube = 37000;   % unit: kPa 
Ep = 33000000;      % unit: kPa 

% settlement parameters
Ms = 0.001;
Ke = 0.4;
Lo = 0;         % unit: m
slim_single = 5;      % unit: mm
slim_group = 10;      % unit: mm
slim_percent_single = 1;   % settlement for single pile limit (% * pile diameter)
slim_percent_group = 2;   % settlement for group pile limit (% * pile diameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input tables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pile and basement
pile_dim = {B, Lmin, Lmax, Lint}; 
basement_data = [column_width_span, basement_depth, column_length_span, soil_otherside_depth];
ciria_data = {Mmin, Mmax, Space};
%% STR/GEO
A1_STR_GEO = [A1_STR_GEO_Dead_Unf, A1_STR_GEO_Dead_F, A1_STR_GEO_Live, A1_STR_GEO_wind]; % load STR/GEO Design combination 1
A2_STR_GEO = [A2_STR_GEO_Dead_Unf, A2_STR_GEO_Dead_F, A2_STR_GEO_Live, A2_STR_GEO_wind]; % load STR/GEO Design combination 2
M1_STR_GEO = [M1_STR_GEO_phi, M1_STR_GEO_c, M1_STR_GEO_Cu, M1_STR_GEO_den]; % material STR/GEO Design combination 1
M2_STR_GEO = [M2_STR_GEO_phi, M2_STR_GEO_c, M2_STR_GEO_Cu, M2_STR_GEO_den]; % material STR/GEO Design combination 2
R1_STR_GEO = [R1_STR_GEO_base, R1_STR_GEO_sh_comp, R1_STR_GEO_sh_ten, R1_STR_GEO_comb_comp];  % resistance STR/GEO Design combination 1
R4_STR_GEO = [R4_STR_GEO_base, R4_STR_GEO_sh_comp, R4_STR_GEO_sh_ten, R4_STR_GEO_comb_comp];  % resistance STR/GEO Design combination 2
STR_GEO_AMR = {A1_STR_GEO, A2_STR_GEO, M1_STR_GEO, M2_STR_GEO, R1_STR_GEO, R4_STR_GEO};
%% design parameters 
Beta_data = [Ks, del_fill, del_tergravel];
conc_data = [fck_cube, Ep];
sett_data = [Ms, Ke, Lo, slim_single, slim_group, slim_percent_single, slim_percent_group];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,v]=size(column_list);
for sc = 1:v
    
    %% Inputs
    column_location = column_list{sc};
    P_list = [Dead_P_list(sc), Live_P_list(sc), Wind_P_list(sc)];
    soil = dlmread(soil_list{sc},',',0,0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STR/GEO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sorted_ULSnSLS_Design, ULSnSLS_Design] = Z5_ULSnSLS_CFA_wall(P_list, ciria_data, STR_GEO_AMR, Beta_data, conc_data, soil, basement_data, pile_dim, sett_data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% design output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ULS and SLS - Z5
    % dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_C1_comp_',column_location,'.csv'),ULSnSLS_Design{1},'delimiter',',')
    % dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_C2_comp_',column_location,'.csv'),ULSnSLS_Design{2},'delimiter',',')
    % dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_C1_ten_',column_location,'.csv'),ULSnSLS_Design{3},'delimiter',',')
    % dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_C2_ten_',column_location,'.csv'),ULSnSLS_Design{4},'delimiter',',')
    dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_overall_comp_',column_location,'.csv'),ULSnSLS_Design{5},'delimiter',',')
    dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_overall_ten_',column_location,'.csv'),ULSnSLS_Design{6},'delimiter',',')
    dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_V sorted_comp_',column_location,'.csv'),sorted_ULSnSLS_Design{1},'delimiter',',')
    dlmwrite(strcat('ZZ5_ULSnSLS CFA wall design_V sorted_ten_',column_location,'.csv'),sorted_ULSnSLS_Design{2},'delimiter',',')
    
end

disp('done')
