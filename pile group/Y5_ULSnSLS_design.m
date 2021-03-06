% Copyright (C) 2016 Enok Cheon

% Author: Enok Cheon <Enok Cheon@LAPTOP-I5VE6B9I>
% Created: 2016-12-30

function [ULSnSLS_Group_Design, every_results, number_of_pile_types_ULSnSLS_Group_Design, ULSnSLS_Group_Cost] = Y5_ULSnSLS_design(P_list, Mx_list, My_list, soil, pile_dim, basement_data, Beta_data, conc_data, Nlim, sett_data, STR_GEO_AMR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% design parameters
% loads
% axial load
Dead_P = P_list(1);
Live_P = P_list(2);
Wind_P = P_list(3);
% moment x-x direction load
Dead_Mx = Mx_list(1);
Live_Mx = Mx_list(2);
Wind_Mx = Mx_list(3);
% moment y-y direction load
Dead_My = My_list(1);
Live_My = My_list(2);
Wind_My = My_list(3);

% pile dimension limits
D = pile_dim{1};
Lmin = pile_dim{2};
Lmax = pile_dim{3};
Lint = pile_dim{4};

% basement data
basement_depth = basement_data(2);

%% analysis parameters
% beta method
Ks = Beta_data(1);
del_fill = Beta_data(2);
del_tergravel = Beta_data(3);

% settlement data
Ms = sett_data(1);
Ke = sett_data(2);
Lo = sett_data(3);
slim_single = sett_data(4);
slim_group = sett_data(5);
slim_percent_single = sett_data(6);
slim_percent_group = sett_data(7); 
Nlim = sett_data(8);

% concrete data
fck_cube = conc_data(1);
Ep = conc_data(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLS calculation iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
comb = {[1, 3, 5],[2, 4, 6]};		% factors to extract from STR_GEO_AMR. (1) is Comb 1 and (2) is Comb 2

% matrix to store analysed values
ULSnSLS_design_C1_comp = [];
ULSnSLS_design_C2_comp = [];
ULSnSLS_design_C1_ten = [];
ULSnSLS_design_C2_ten = [];
ULSnSLS_Design_comp = [];
ULSnSLS_Design_ten = [];
sorted_ULSnSLS_Design_comp = [];
sorted_ULSnSLS_Design_ten = [];
tab = 1;
every_results = [];

for i = 1:2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% soil
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% [depth, length, density, undrained shear strength, friction angle, cohesion, alpha, E, v (poisson ratio)]
	ULSnSLS_soil = [soil(:,1), soil(:,2), soil(:,3)/STR_GEO_AMR{comb{i}(2)}(4), soil(:,4)/STR_GEO_AMR{comb{i}(2)}(3), ...
	            atand(tand(soil(:,5))/STR_GEO_AMR{comb{i}(2)}(1)), soil(:,6)/STR_GEO_AMR{comb{i}(2)}(2), soil(:,7), soil(:,8), soil(:,9)];
	ULSnSLS_soil_fill = ULSnSLS_soil(1,:);
	ULSnSLS_soil_orgclay = ULSnSLS_soil(2,:);
	ULSnSLS_soil_tergravel = ULSnSLS_soil(3,:);
	ULSnSLS_soil_Lonclay1 = ULSnSLS_soil(4,:);
	ULSnSLS_soil_Lonclay2 = ULSnSLS_soil(5,:);
	ULSnSLS_soil_Lonclay3 = ULSnSLS_soil(6,:);

	row = 1;
    for d = D
        A = 0.25*pi*d.^2;         % shaft cross-sectional area
        for Z = [Lmin+basement_depth: Lint: Lmax+basement_depth]
            L = Z-basement_depth;
            for N = 3:Nlim
        	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		%% input from the ULS analysis - N 
        		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % assumption: centre-to-centre distance between piles (dx) is 3 x diameter
                switch N
                    case 3  % equilateral triangular arrangement
                        dx = 3*d;                   % space between piles
                        % distribution of moment loads
                        xn = [-dx/2; 0; -dx/2]; 
                        yn = [-0.25*sqrt(3)*dx; 0.25*sqrt(3)*dx; -0.25*sqrt(3)*dx];                
                        Ixx = (9/16)*dx^2;
                        Iyy = 0.5*dx^2;
                        % group
                        bb = 0.5*sqrt(3)*(dx+d);    % width (base)
                        ll = dx+d;                  % length (base)
                        Ab = bb*ll*0.5;             % base area
                        perimeter = 2*(d+2*dx);     
                        As = L*perimeter;           % shaft area    
                        % distribution of axial loads
                        N_pile_I = 0;               % number of piles interior - surrounded by piles
                        N_pile_E = 0;               % number of piles on exterior that is not corner
                        N_pile_C = 3;               % number of piles at corner
                        pile_coordinates = [0,0; dx/2,bb; dx,0];        % (x,y) coordinate showing pile locations - left->right, down->up
                        pile_types = [3; 3; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C
                    case 4  % 2x2 layout
                        dx = 3*d;
                        % distribution of axial loads
                        N_pile_I = 0;
                        N_pile_E = 0;
                        N_pile_C = 4;
                        pile_coordinates = [0,0; 0,dx; dx,dx; dx,0];  
                        pile_types = [3; 3; 3; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C  
                        % distribution of moment loads
                        xn = [-dx/2; -dx/2; dx/2; dx/2];  
                        yn = [-dx/2; dx/2; dx/2; -dx/2];   
                        Ixx = dx^2;
                        Iyy = dx^2;
                        % group
                        bb = dx+d;
                        ll = dx+d;
                        Ab = (d+dx)^2;       
                        perimeter = 4*(d+dx); 
                        As = L*perimeter;
                    case 5  % (2-1-2) layout rectangular with a pile in centre
                        dx = 3*d;
                        % distribution of moment loads
                        xn = [-dx; -dx; 0; dx; dx];  
                        yn = [-dx; dx; 0; dx; -dx];  
                        Ixx = 4*dx^2;
                        Iyy = 4*dx^2;
                        % group
                        bb = (2*dx+d);
                        ll = (2*dx+d);
                        Ab = (2*dx+d)^2;     
                        perimeter = 4*(2*dx+d); 
                        As = L*perimeter;
                        % distribution of axial loads
                        N_pile_I = 1;
                        N_pile_E = 0;
                        N_pile_C = 4;
                        pile_coordinates = [0,0; 2*dx,0; dx,dx; 2*dx,0; 2*dx,2*dx];
                        pile_types = [3; 3; 1; 3; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C        
                    case 6   % 2x3 layout 
                        dx = 3*d;
                        % distribution of moment loads
                        xn = [-dx; -dx; 0; 0; dx; dx];
                        yn = [-dx/2; dx/2; -dx/2; dx/2; -dx/2; dx/2];
                        Ixx = 4*dx^2;
                        Iyy = 1.5*dx^2;
                        % group
                        bb = dx+d;
                        ll = (d+2*dx);
                        Ab = bb*ll;      
                        perimeter = 2*((d+dx)+(d+2*dx)); 
                        As = L*perimeter;
                        % distribution of axial loads
                        N_pile_I = 0;
                        N_pile_E = 2;
                        N_pile_C = 4;
                        pile_coordinates = [0,0; 0,dx; dx,0; dx,dx; 2*dx,0; 2*dx,dx];        % coordinate showing pile locations
                        pile_types = [3; 3; 2; 2; 3; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C
                    case 7  % (2-3-2) layout, regular hexagon with centre
                        dx = 3*d;
                        % distribution of moment loads
                        xn = [-dx; -dx/2; -dx/2; 0; dx/2; dx/2; dx];
                        yn = [0; -0.5*sqrt(3)*dx; 0.5*sqrt(3)*dx; 0; -0.5*sqrt(3)*dx; 0.5*sqrt(3)*dx; 0];
                        Ixx = 3*dx^2;
                        Iyy = 3*dx^2;
                        % group
                        bb = sqrt(3)*(dx+d);
                        ll = (d+2*dx);
                        Ab = 1.5*sqrt(3)*(dx+d)^2;           
                        perimeter = 6*(dx+d); 
                        As = L*perimeter;
                        % distribution of axial loads
                        N_pile_I = 1;
                        N_pile_E = 4;
                        N_pile_C = 2;
                        pile_coordinates = [0,0; dx/2,-0.5*sqrt(3)*dx; dx/2,0.5*sqrt(3)*dx; dx,0; dx*1.5,-0.5*sqrt(3)*dx; dx*1.5,0.5*sqrt(3)*dx; 2*dx,0];        % coordinate showing pile locations
                        pile_types = [3; 2; 2; 1; 2; 2; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C
                    case 8  % 2x4 layout
                        dx = 3*d;
                        % distribution of moment loads
                        xn = [-1.5*dx; -1.5*dx; -0.5*dx; -0.5*dx; 0.5*dx; 0.5*dx; 1.5*dx; 1.5*dx];
                        yn = [-dx/2; dx/2; -dx/2; dx/2; -dx/2; dx/2; -dx/2; dx/2];
                        Ixx = 2*dx^2;
                        Iyy = 10*dx^2;
                        % group
                        bb = dx+d;
                        ll = (d+3*dx);
                        Ab = bb*ll;
                        perimeter = 2*(bb+ll); 
                        As = L*perimeter;
                        % distribution of axial loads
                        N_pile_I = 0;
                        N_pile_E = 4;
                        N_pile_C = 4;
                        pile_coordinates = [0,0; 0,dx; dx,0; dx,dx; 2*dx,0; 2*dx,dx; 3*dx,0; 3*dx,dx];      % coordinate showing pile locations
                        pile_types = [3; 3; 2; 2; 2; 2; 3; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C
                    case 9  % 3x3 layout
                        dx = 3*d;
                        % distribution of moment loads
                        xn = [-dx; -dx; -dx; 0; 0; 0; dx; dx; dx];
                        yn = [-dx; 0; dx; -dx; 0; dx; -dx; 0; dx];
                        Ixx = 6*dx^2;
                        Iyy = 6*dx^2;
                        % group
                        bb = (d+2*dx);
                        ll = (d+2*dx);
                        Ab = (d+2*dx)^2;
                        perimeter = 4*(d+2*dx); 
                        As = L*perimeter;
                        % distribution of axial loads
                        N_pile_I = 1;
                        N_pile_E = 4;
                        N_pile_C = 4;
                        pile_coordinates = [0,0; 0,dx; 0,2*dx; dx,0; dx,dx; dx,2*dx; 2*dx,0; 2*dx,dx; 2*dx,2*dx];        % coordinate showing pile locations
                        pile_types = [3; 2; 3; 2; 1; 2; 3; 2; 3];     % 1 = N_pile_I, 2 = N_pile_E, 3 = N_pile_C
                    otherwise
                        error('too many piles or not enough - check your codes')
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% load combination analysis - total axial loads
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% factored loads
                ULS_Dead_Unf = STR_GEO_AMR{comb{i}(1)}(1)*[Dead_P, Dead_Mx, Dead_My];
                ULS_Dead_F = STR_GEO_AMR{comb{i}(1)}(2)*[Dead_P, Dead_Mx, Dead_My];
                ULS_Live_Unf = STR_GEO_AMR{comb{i}(1)}(3)*[Live_P, Live_Mx, Live_My];
                ULS_Wind_Unf = STR_GEO_AMR{comb{i}(1)}(4)*[Wind_P, Wind_Mx, Wind_My];
                SLS_Dead_Unf = [Dead_P, Dead_Mx, Dead_My];
                SLS_Dead_F = [Dead_P, Dead_Mx, Dead_My];
                SLS_Live_Unf = [Live_P, Live_Mx, Live_My];
                SLS_Wind_Unf = [Wind_P, Wind_Mx, Wind_My];

                ULS_Load_comp_total = 0;
                ULS_Load_ten_total = 0;

                %% dead load only
                % combination case 1
                load_combination_case = [];
                if Dead_P == 0 
                    error('you seem to have no dead load from your structures')
                elseif Dead_P ~= 0                        % construction load
                    if Dead_P > 0
                        ULS_load_1_comp_total = ULS_Dead_Unf(1);   % for overall group pile load
                        ULS_load_1_ten_total = 0;
                        load_combination_case(1) = 19;
                    else 
                        ULS_load_1_comp_total =  0; % for overall group pile load
                        ULS_load_1_ten_total = ULS_Dead_Unf(1);
                        load_combination_case(1) = 91;
                    end
                    ULS_M1 = [ULS_Dead_Unf(2), ULS_Dead_Unf(3)]; 
                end
                %% dead and live load (no wind)
                if (Live_P == 0 && Dead_P > 0)                                       % no live
                    ULS_load_2_comp_total = 0;
                    ULS_load_2_ten_total = 0;
                    ULS_M2 = [0, 0];
                    load_combination = 0;
                % combination case 2
                elseif (Live_P > 0 && Dead_P > 0) || (Live_P< 0 && Dead_P < 0)                                  % gravity only - both compression/tension
                    ULS_load_2 = ULS_Dead_Unf(1) + ULS_Live_Unf(1);
                    if ULS_load_2 > 0
                        ULS_load_2_comp_total = ULS_load_2;
                        ULS_load_2_ten_total = 0;
                        load_combination = 29;
                    else 
                        ULS_load_2_comp_total = 0;
                        ULS_load_2_ten_total = ULS_load_2;
                        load_combination = 92;
                    end
                    ULS_M2 = [ULS_Dead_Unf(2)+ ULS_Live_Unf(2),ULS_Dead_Unf(3)+ULS_Live_Unf(3)]; 
                % combination case 1, 3
                elseif (Live_P  > 0 && Dead_P < 0) || (Live_P< 0 && Dead_P > 0)         % gravity only - different load type for each loads
                    % axial
                    ULS_load_2a = ULS_Dead_Unf(1);                          % live is favourable; therefore, should not be considered
                    ULS_load_2b = ULS_Dead_F(1) + ULS_Live_Unf(1);          % dead is favourable
                    ULS_load_2_comp_total = max(ULS_load_2a,ULS_load_2b);
                    ULS_load_2_ten_total = min(ULS_load_2a,ULS_load_2b);
                    % moment
                    ULS_M_2a = [ULS_Dead_Unf(2), ULS_Dead_Unf(3)];
                    ULS_M_2b = [ULS_Dead_F(2) + ULS_Live_Unf(2), ULS_Dead_F(3) + ULS_Live_Unf(3)];
                    ULS_M2 = [];
                    if abs(ULS_M_2a(1)) >= abs(ULS_M_2b(1))
                        ULS_M2(1) = ULS_M_2a(1);
                    else
                        ULS_M2(1) = ULS_M_2b(1);
                    end
                    if abs(ULS_M_2a(2)) >= abs(ULS_M_2b(2))
                        ULS_M2(2) = ULS_M_2a(2);
                    else
                        ULS_M2(2) = ULS_M_2b(2);
                    end
                    % orientation
                    if ULS_load_2_comp_total == ULS_load_2a && ULS_load_2_ten_total == ULS_load_2b 
                        load_combination = 13;
                    elseif ULS_load_2_comp_total == ULS_load_2b && ULS_load_2_ten_total == ULS_load_2a 
                        load_combination = 31;
                    end
                end
                load_combination_case(2) = load_combination;
                %% all load combinations
                if Wind_P == 0 && Dead_P > 0                    % no wind
                    ULS_load_3_comp_total = 0;
                    ULS_load_3_ten_total = 0;
                    ULS_M3 = [0, 0];
                    load_combination = 0;
                % combination case 4
                elseif (Live_P> 0 && Dead_P > 0 && Wind_P > 0) || (Live_P< 0 && Dead_P < 0 && Wind_P < 0)  % all in compression/tension
                    ULS_load_3 =  ULS_Dead_Unf(1) + ULS_Live_Unf(1) + ULS_Wind_Unf(1);
                    if ULS_load_3 > 0
                        ULS_load_3_comp_total = ULS_load_3;
                        ULS_load_3_ten_total = 0;
                        load_combination = 49;
                    else 
                        ULS_load_3_comp_total = 0;
                        ULS_load_3_ten_total = ULS_load_3;
                        load_combination = 94;
                    end            
                    ULS_M3 = [ULS_Dead_Unf(2)+ULS_Live_Unf(2)+ULS_Wind_Unf(2), ULS_Dead_Unf(3)+ULS_Live_Unf(3)+ULS_Wind_Unf(3)]; 
                % combination case 5,1
                elseif (Live_P> 0 && Dead_P < 0 && Wind_P > 0)||(Live_P< 0 && Dead_P > 0 && Wind_P < 0)   % live and wind are same; dead is different 
                    ULS_load_3a = ULS_Dead_F(1) + ULS_Live_Unf(1) + ULS_Wind_Unf(1);
                    ULS_load_3b = ULS_Dead_Unf(1);
                    ULS_load_3_comp_total = max(ULS_load_3a,ULS_load_3b);
                    ULS_load_3_ten_total = min(ULS_load_3a,ULS_load_3b);

                    ULS_M_3a = [ULS_Dead_F(2) + ULS_Live_Unf(2) + ULS_Wind_Unf(2), ULS_Dead_F(3) + ULS_Live_Unf(3) + ULS_Wind_Unf(3)];
                    ULS_M_3b = [ULS_Dead_Unf(2), ULS_Dead_Unf(3)];
                    ULS_M3 = [];
                    if abs(ULS_M_3a(1)) >= abs(ULS_M_3b(1))
                        ULS_M3(1) = ULS_M_3a(1);
                    else
                        ULS_M3(1) = ULS_M_3b(1);
                    end
                    if abs(ULS_M_3a(2)) >= abs(ULS_M_3b(2))
                        ULS_M3(2) = ULS_M_3a(2);
                    else
                        ULS_M3(2) = ULS_M_3b(2);
                    end

                    if ULS_load_3_comp_total == ULS_load_3a && ULS_load_3_ten_total == ULS_load_3b 
                        load_combination = 51;
                    elseif ULS_load_3_comp_total == ULS_load_3b && ULS_load_3_ten_total == ULS_load_3a 
                        load_combination = 15;
                    end
                % combination case 2
                elseif (Live_P> 0 && Dead_P > 0 && Wind_P < 0)||(Live_P < 0 && Dead_P < 0 && Wind_P > 0)   % live and dead are same; wind is different
                    ULS_load_3 = ULS_Dead_Unf(1) + ULS_Live_Unf(1);
                    if ULS_load_3 > 0
                        ULS_load_3_comp_total = ULS_load_3;
                        ULS_load_3_ten_total = 0;
                        load_combination = 29;
                    else 
                        ULS_load_3_comp_total = 0;
                        ULS_load_3_ten_total = ULS_load_3;
                        load_combination = 92;
                    end
                    ULS_M3 = [ULS_Dead_Unf(2)+ULS_Live_Unf(2), ULS_Dead_Unf(3)+ULS_Live_Unf(3)];
                % combination case 6,3
                elseif (Live_P < 0 && Dead_P > 0 && Wind_P > 0)||(Live_P> 0 && Dead_P < 0 && Wind_P < 0)   % dead and wind are same; live is different
                    ULS_load_3a = ULS_Dead_Unf(1) + ULS_Wind_Unf(1);
                    ULS_load_3b = ULS_Dead_F(1) + ULS_Live_Unf(1);  
                    ULS_load_3_comp_total = max(ULS_load_3a,ULS_load_3b);
                    ULS_load_3_ten_total = min(ULS_load_3a,ULS_load_3b);

                    ULS_M_3a = [ULS_Dead_Unf(2) + ULS_Wind_Unf(2), ULS_Dead_Unf(3) + ULS_Wind_Unf(3)];
                    ULS_M_3b = [ULS_Dead_F(2) + ULS_Live_Unf(2), ULS_Dead_F(3) + ULS_Live_Unf(3)];
                    ULS_M3 = [];
                    if abs(ULS_M_3a(1)) >= abs(ULS_M_3b(1))
                        ULS_M3(1) = ULS_M_3a(1);
                    else
                        ULS_M3(1) = ULS_M_3b(1);
                    end
                    if abs(ULS_M_3a(2)) >= abs(ULS_M_3b(2))
                        ULS_M3(2) = ULS_M_3a(2);
                    else
                        ULS_M3(2) = ULS_M_3b(2);
                    end

                    if ULS_load_3_comp_total == ULS_load_3a && ULS_load_3_ten_total == ULS_load_3b 
                        load_combination = 63;
                    elseif ULS_load_3_comp_total == ULS_load_3b && ULS_load_3_ten_total == ULS_load_3a 
                        load_combination = 36;
                    end
                end
                load_combination_case(3) = load_combination;

                % compile axial loads to design 
                ULS_Load_comp_total = max([ULS_load_1_comp_total, ULS_load_2_comp_total, ULS_load_3_comp_total]);
                ULS_Load_ten_total = min([ULS_load_1_ten_total, ULS_load_2_ten_total, ULS_load_3_ten_total]);
                % compile moment loads to design
                % Mx
                if abs(ULS_M1(1)) == max([abs(ULS_M1(1)), abs(ULS_M2(1)), abs(ULS_M3(1))])
                    ULS_Load_Mx = ULS_M1(1);
                elseif abs(ULS_M2(1)) == max([abs(ULS_M1(1)), abs(ULS_M2(1)), abs(ULS_M3(1))])
                    ULS_Load_Mx = ULS_M2(1);
                elseif abs(ULS_M3(1)) == max([abs(ULS_M1(1)), abs(ULS_M2(1)), abs(ULS_M3(1))])
                    ULS_Load_Mx = ULS_M3(1);
                end
                % My
                if abs(ULS_M1(2)) == max([abs(ULS_M1(2)), abs(ULS_M2(2)), abs(ULS_M3(2))])
                    ULS_Load_My = ULS_M1(2);
                elseif abs(ULS_M2(2)) == max([abs(ULS_M1(2)), abs(ULS_M2(2)), abs(ULS_M3(2))])
                    ULS_Load_My = ULS_M2(2);
                elseif abs(ULS_M3(2)) == max([abs(ULS_M1(2)), abs(ULS_M2(2)), abs(ULS_M3(2))])
                    ULS_Load_My = ULS_M3(2);
                end

                % determine the critical load combination
                if (ULS_Load_comp_total ~= 0 && ULS_Load_ten_total == 0)
                    if ULS_Load_comp_total == ULS_load_1_comp_total 
                        critical_load_combi = load_combination_case(1);
                    elseif ULS_Load_comp_total == ULS_load_2_comp_total
                        critical_load_combi = load_combination_case(2);
                    elseif ULS_Load_comp_total == ULS_load_3_comp_total
                        critical_load_combi = load_combination_case(3);
                    end
                elseif (ULS_Load_comp_total == 0 && ULS_Load_ten_total ~= 0)
                    if ULS_Load_ten_total == ULS_load_1_ten_total 
                        critical_load_combi = load_combination_case(1);
                    elseif ULS_Load_ten_total == ULS_load_2_ten_total
                        critical_load_combi = load_combination_case(2);
                    elseif ULS_Load_ten_total == ULS_load_3_ten_total
                        critical_load_combi = load_combination_case(3);
                    end
                elseif ULS_Load_comp_total ~= 0 && ULS_Load_ten_total ~= 0
                    if ULS_Load_ten_total == ULS_load_1_ten_total && ULS_Load_comp_total == ULS_load_1_comp_total 
                        critical_load_combi = load_combination_case(1);
                    elseif ULS_Load_ten_total == ULS_load_2_ten_total && ULS_Load_comp_total == ULS_load_2_comp_total
                        critical_load_combi = load_combination_case(2);
                    elseif ULS_Load_ten_total == ULS_load_3_ten_total && ULS_Load_comp_total == ULS_load_3_comp_total
                        critical_load_combi = load_combination_case(3);
                    end
                end

                % list of total axial loads to design
                if isempty(ULS_Load_comp_total) == 1 || isempty(ULS_Load_ten_total) == 1 % no loads
                    error('something wrong with the loads')
                elseif (isempty(ULS_Load_comp_total) == 0 && isempty(ULS_Load_ten_total) == 0)  % both tension and compression occurs
                    ULS_Load_total_list = [ULS_Load_comp_total, ULS_Load_ten_total];
                end

                %% SLS load -allowable loads
                switch critical_load_combi
                    case 19
                        SLS_Load_comp_total = SLS_Dead_Unf(1);
                        SLS_Load_ten_total = 0;
                        SLS_Load_Mx = SLS_Dead_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3);
                    case 91
                        SLS_Load_comp_total = 0;
                        SLS_Load_ten_total = SLS_Dead_Unf(1); 
                        SLS_Load_Mx = SLS_Dead_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3);
                    case 13
                        SLS_Load_comp_total = SLS_Dead_Unf(1); 
                        SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
                        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
                    case 31
                        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
                        SLS_Load_ten_total = SLS_Dead_Unf(1); 
                        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
                    case 15
                        SLS_Load_comp_total = SLS_Dead_Unf(1); 
                        SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
                        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3));
                    case 51
                        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
                        SLS_Load_ten_total = SLS_Dead_Unf(1); 
                        SLS_Load_Mx = max(SLS_Dead_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3));
                    case 29
                        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
                        SLS_Load_ten_total = 0;
                        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3);
                    case 92
                        SLS_Load_comp_total = 0;
                        SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
                        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3);
                    case 49
                        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
                        SLS_Load_ten_total = 0;
                        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3);
                    case 94
                        SLS_Load_comp_total = 0;
                        SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
                        SLS_Load_Mx = SLS_Dead_Unf(2)+SLS_Live_Unf(2)+SLS_Wind_Unf(2);
                        SLS_Load_My = SLS_Dead_Unf(3)+SLS_Live_Unf(3)+SLS_Wind_Unf(3);
                    case 36
                        SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
                        SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
                        SLS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
                    case 63
                        SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
                        SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
                        SLS_Load_Mx = max(SLS_Dead_Unf(2)+SLS_Wind_Unf(2), SLS_Dead_F(2)+SLS_Live_Unf(2));
                        SLS_Load_My = max(SLS_Dead_Unf(3)+SLS_Wind_Unf(3), SLS_Dead_F(3)+SLS_Live_Unf(3));
                end

                % list of total axial loads to design
                if isempty(SLS_Load_comp_total) == 1 || isempty(SLS_Load_ten_total) == 1 % no loads
                    error('something wrong with the loads')
                elseif (isempty(SLS_Load_comp_total) == 0 && isempty(SLS_Load_ten_total) == 0)  % both tension and compression occurs
                    SLS_Load_total_list = [SLS_Load_comp_total, SLS_Load_ten_total];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% iterate ULS_Group_Design
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for loads_cases = 1:2
                    ULS_load_total = ULS_Load_total_list(loads_cases);
                    SLS_load_total = SLS_Load_total_list(loads_cases);
                    if ULS_load_total == 0 || SLS_load_total == 0
                        if loads_cases == 1
                            continue
                        elseif loads_cases == 2
                            break
                        end
                    end

            		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% ULS additional load on each piles based on moments
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                    ULS_loads_from_M = []; 
                    for m = 1:length(pile_types)
                        if abs(ULS_Load_Mx) >= abs(ULS_Load_My)
                            if Ixx >= Iyy
                                ULS_loads_from_M1 = (ULS_Load_Mx*yn(m))/Ixx + (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M2 = (ULS_Load_Mx*yn(m))/Ixx - (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M3 = -(ULS_Load_Mx*yn(m))/Ixx + (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M4 = -(ULS_Load_Mx*yn(m))/Ixx - (ULS_Load_My*xn(m))/Iyy;
                                orientation = 0;
                            elseif Ixx < Iyy
                                ULS_loads_from_M1 = (ULS_Load_Mx*yn(m))/Iyy + (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M2 = (ULS_Load_Mx*yn(m))/Iyy - (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M3 = -(ULS_Load_Mx*yn(m))/Iyy + (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M4 = -(ULS_Load_Mx*yn(m))/Iyy - (ULS_Load_My*xn(m))/Ixx;
                                orientation = 1;
                            end
                        elseif abs(ULS_Load_Mx) < abs(ULS_Load_My)
                            if Ixx > Iyy
                                ULS_loads_from_M1 = (ULS_Load_Mx*yn(m))/Ixx + (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M2 = (ULS_Load_Mx*yn(m))/Ixx - (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M3 = -(ULS_Load_Mx*yn(m))/Ixx + (ULS_Load_My*xn(m))/Iyy;
                                ULS_loads_from_M4 = -(ULS_Load_Mx*yn(m))/Ixx - (ULS_Load_My*xn(m))/Iyy;
                                orientation = 0;
                            elseif Ixx < Iyy
                                ULS_loads_from_M1 = (ULS_Load_Mx*yn(m))/Iyy + (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M2 = (ULS_Load_Mx*yn(m))/Iyy - (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M3 = -(ULS_Load_Mx*yn(m))/Iyy + (ULS_Load_My*xn(m))/Ixx;
                                ULS_loads_from_M4 = -(ULS_Load_Mx*yn(m))/Iyy - (ULS_Load_My*xn(m))/Ixx;
                                orientation = 1;
                            end
                        end
                                                 
                        if ULS_Load_Mx == 0 && ULS_Load_My == 0
                            ULS_loads_from_M(m) = 0;
                        elseif ULS_load_total > 0 && (ULS_Load_Mx ~= 0 || ULS_Load_My ~= 0)
                            if ULS_loads_from_M1 == max([ULS_loads_from_M1,ULS_loads_from_M2,ULS_loads_from_M3,ULS_loads_from_M4]) && ULS_loads_from_M1 > 0
                                ULS_loads_from_M(m) = ULS_loads_from_M1;
                            elseif ULS_loads_from_M2 == max([ULS_loads_from_M1,ULS_loads_from_M2,ULS_loads_from_M3,ULS_loads_from_M4]) && ULS_loads_from_M2 > 0
                                ULS_loads_from_M(m) = ULS_loads_from_M2;
                            elseif ULS_loads_from_M3 == max([ULS_loads_from_M1,ULS_loads_from_M2,ULS_loads_from_M3,ULS_loads_from_M4]) && ULS_loads_from_M3 > 0
                                ULS_loads_from_M(m) = ULS_loads_from_M3;
                            elseif ULS_loads_from_M4 == max([ULS_loads_from_M1,ULS_loads_from_M2,ULS_loads_from_M3,ULS_loads_from_M4]) && ULS_loads_from_M4 > 0
                                ULS_loads_from_M(m) = ULS_loads_from_M4;
                            end
                        elseif ULS_load_total < 0 && (ULS_Load_Mx ~= 0 || ULS_Load_My ~= 0)
                            if abs(ULS_loads_from_M1) == min([abs(ULS_loads_from_M1),abs(ULS_loads_from_M2),abs(ULS_loads_from_M3),abs(ULS_loads_from_M4)]) && ULS_loads_from_M1 < 0
                                ULS_loads_from_M(m) = ULS_loads_from_M1;
                            elseif abs(ULS_loads_from_M2) == min([abs(ULS_loads_from_M1),abs(ULS_loads_from_M2),abs(ULS_loads_from_M3),abs(ULS_loads_from_M4)]) && ULS_loads_from_M2 < 0
                                ULS_loads_from_M(m) = ULS_loads_from_M2;
                            elseif abs(ULS_loads_from_M3) == min([abs(ULS_loads_from_M1),abs(ULS_loads_from_M2),abs(ULS_loads_from_M3),abs(ULS_loads_from_M4)]) && ULS_loads_from_M3 < 0
                                ULS_loads_from_M(m) = ULS_loads_from_M3;
                            elseif abs(ULS_loads_from_M4) == min([abs(ULS_loads_from_M1),abs(ULS_loads_from_M2),abs(ULS_loads_from_M3),abs(ULS_loads_from_M4)]) && ULS_loads_from_M4 < 0
                                ULS_loads_from_M(m) = ULS_loads_from_M4;
                            end
                        end
                    end

                    % sort data for each pile type
                    ULS_loads_from_M_per_type = {[], [], []};
                    for m = 1:length(pile_types)
                        if ULS_loads_from_M(pile_types(m)) == 0
                            ULS_loads_from_M_per_type{pile_types(m)} = [ULS_loads_from_M_per_type{pile_types(m)}, 0];
                        else
                            ULS_loads_from_M_per_type{pile_types(m)} = [ULS_loads_from_M_per_type{pile_types(m)}, ULS_loads_from_M(pile_types(m))];
                        end
                    end
                    for m = 1:3
                        if isempty(ULS_loads_from_M_per_type{m}) == 1
                            ULS_loads_from_M_per_type{m} = [0];
                        end
                    end

                    % sort maximum additional load for each pile type
                    if ULS_load_total > 0
                        ULS_load_P_from_M_per_pile_type_I = max(ULS_loads_from_M_per_type{1});
                        ULS_load_P_from_M_per_pile_type_E = max(ULS_loads_from_M_per_type{2});
                        ULS_load_P_from_M_per_pile_type_C = max(ULS_loads_from_M_per_type{3});
                    elseif ULS_load_total < 0 
                        ULS_load_P_from_M_per_pile_type_I = min(ULS_loads_from_M_per_type{1});
                        ULS_load_P_from_M_per_pile_type_E = min(ULS_loads_from_M_per_type{2});
                        ULS_load_P_from_M_per_pile_type_C = min(ULS_loads_from_M_per_type{3});
                    end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% SLS additional load on each piles based on moments
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                    SLS_loads_from_M = []; 
                    for m = 1:length(pile_types)
                        if abs(SLS_Load_Mx) >= abs(SLS_Load_My)
                            if Ixx >= Iyy
                                SLS_loads_from_M1 = (SLS_Load_Mx*yn(m))/Ixx + (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M2 = (SLS_Load_Mx*yn(m))/Ixx - (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M3 = -(SLS_Load_Mx*yn(m))/Ixx + (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M4 = -(SLS_Load_Mx*yn(m))/Ixx - (SLS_Load_My*xn(m))/Iyy;
                                orientation = 0;
                            elseif Ixx < Iyy
                                SLS_loads_from_M1 = (SLS_Load_Mx*yn(m))/Iyy + (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M2 = (SLS_Load_Mx*yn(m))/Iyy - (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M3 = -(SLS_Load_Mx*yn(m))/Iyy + (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M4 = -(SLS_Load_Mx*yn(m))/Iyy - (SLS_Load_My*xn(m))/Ixx;
                                orientation = 1;
                            end
                        elseif abs(SLS_Load_Mx) < abs(SLS_Load_My)
                            if Ixx > Iyy
                                SLS_loads_from_M1 = (SLS_Load_Mx*yn(m))/Ixx + (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M2 = (SLS_Load_Mx*yn(m))/Ixx - (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M3 = -(SLS_Load_Mx*yn(m))/Ixx + (SLS_Load_My*xn(m))/Iyy;
                                SLS_loads_from_M4 = -(SLS_Load_Mx*yn(m))/Ixx - (SLS_Load_My*xn(m))/Iyy;
                                orientation = 0;
                            elseif Ixx < Iyy
                                SLS_loads_from_M1 = (SLS_Load_Mx*yn(m))/Iyy + (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M2 = (SLS_Load_Mx*yn(m))/Iyy - (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M3 = -(SLS_Load_Mx*yn(m))/Iyy + (SLS_Load_My*xn(m))/Ixx;
                                SLS_loads_from_M4 = -(SLS_Load_Mx*yn(m))/Iyy - (SLS_Load_My*xn(m))/Ixx;
                                orientation = 1;
                            end
                        end
                                                 
                        if SLS_Load_Mx == 0 && SLS_Load_My == 0
                            SLS_loads_from_M(m) = 0;
                        elseif SLS_load_total > 0 && (SLS_Load_Mx ~= 0 || SLS_Load_My ~= 0)
                            if SLS_loads_from_M1 == max([SLS_loads_from_M1,SLS_loads_from_M2,SLS_loads_from_M3,SLS_loads_from_M4]) && SLS_loads_from_M1 > 0
                                SLS_loads_from_M(m) = SLS_loads_from_M1;
                            elseif SLS_loads_from_M2 == max([SLS_loads_from_M1,SLS_loads_from_M2,SLS_loads_from_M3,SLS_loads_from_M4]) && SLS_loads_from_M2 > 0
                                SLS_loads_from_M(m) = SLS_loads_from_M2;
                            elseif SLS_loads_from_M3 == max([SLS_loads_from_M1,SLS_loads_from_M2,SLS_loads_from_M3,SLS_loads_from_M4]) && SLS_loads_from_M3 > 0
                                SLS_loads_from_M(m) = SLS_loads_from_M3;
                            elseif SLS_loads_from_M4 == max([SLS_loads_from_M1,SLS_loads_from_M2,SLS_loads_from_M3,SLS_loads_from_M4]) && SLS_loads_from_M4 > 0
                                SLS_loads_from_M(m) = SLS_loads_from_M4;
                            end
                        elseif SLS_load_total < 0 && (SLS_Load_Mx ~= 0 || SLS_Load_My ~= 0)
                            if abs(SLS_loads_from_M1) == min([abs(SLS_loads_from_M1),abs(SLS_loads_from_M2),abs(SLS_loads_from_M3),abs(SLS_loads_from_M4)]) && SLS_loads_from_M1 < 0
                                SLS_loads_from_M(m) = SLS_loads_from_M1;
                            elseif abs(SLS_loads_from_M2) == min([abs(SLS_loads_from_M1),abs(SLS_loads_from_M2),abs(SLS_loads_from_M3),abs(SLS_loads_from_M4)]) && SLS_loads_from_M2 < 0
                                SLS_loads_from_M(m) = SLS_loads_from_M2;
                            elseif abs(SLS_loads_from_M3) == min([abs(SLS_loads_from_M1),abs(SLS_loads_from_M2),abs(SLS_loads_from_M3),abs(SLS_loads_from_M4)]) && SLS_loads_from_M3 < 0
                                SLS_loads_from_M(m) = SLS_loads_from_M3;
                            elseif abs(SLS_loads_from_M4) == min([abs(SLS_loads_from_M1),abs(SLS_loads_from_M2),abs(SLS_loads_from_M3),abs(SLS_loads_from_M4)]) && SLS_loads_from_M4 < 0
                                SLS_loads_from_M(m) = SLS_loads_from_M4;
                            end
                        end
                    end

                    % sort data for each pile type
                    SLS_loads_from_M_per_type = {[], [], []};
                    for m = 1:length(pile_types)
                        if SLS_loads_from_M(pile_types(m)) == 0
                            SLS_loads_from_M_per_type{pile_types(m)} = [SLS_loads_from_M_per_type{pile_types(m)}, 0];
                        else
                            SLS_loads_from_M_per_type{pile_types(m)} = [SLS_loads_from_M_per_type{pile_types(m)}, SLS_loads_from_M(pile_types(m))];
                        end
                    end
                    for m = 1:3
                        if isempty(SLS_loads_from_M_per_type{m}) == 1
                            SLS_loads_from_M_per_type{m} = [0];
                        end
                    end

                    % sort maximum additional load for each pile type
                    if SLS_load_total > 0
                        SLS_load_P_from_M_per_pile_type_I = max(SLS_loads_from_M_per_type{1});
                        SLS_load_P_from_M_per_pile_type_E = max(SLS_loads_from_M_per_type{2});
                        SLS_load_P_from_M_per_pile_type_C = max(SLS_loads_from_M_per_type{3});
                    elseif SLS_load_total < 0 
                        SLS_load_P_from_M_per_pile_type_I = min(SLS_loads_from_M_per_type{1});
                        SLS_load_P_from_M_per_pile_type_E = min(SLS_loads_from_M_per_type{2});
                        SLS_load_P_from_M_per_pile_type_C = min(SLS_loads_from_M_per_type{3});
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% ULS - load distribution from group pile
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% soil Young's and Shear Modulus
                    if Z >= ULSnSLS_soil_Lonclay3(1) && Z < ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2) 
                        Eb = ULSnSLS_soil_Lonclay3(8)*1000;
                        v = ULSnSLS_soil_Lonclay3(9);
                        G = Eb/(2*(1+v));
                    elseif Z >= ULSnSLS_soil_Lonclay2(1) && Z < ULSnSLS_soil_Lonclay3(1) 
                        Eb = ULSnSLS_soil_Lonclay2(8)*1000;
                        v = ULSnSLS_soil_Lonclay2(9);
                        G = Eb/(2*(1+v));
                    elseif Z >= ULSnSLS_soil_Lonclay1(1) && Z < ULSnSLS_soil_Lonclay2(1) 
                        Eb = ULSnSLS_soil_Lonclay1(8)*1000;
                        v = ULSnSLS_soil_Lonclay1(9);
                        G = Eb/(2*(1+v));
                    elseif Z >= ULSnSLS_soil_tergravel(1) && Z < ULSnSLS_soil_Lonclay1(1) 
                        Eb = ULSnSLS_soil_tergravel(8)*1000;
                        v = ULSnSLS_soil_tergravel(9);
                        G = Eb/(2*(1+v));
                    elseif Z >= ULSnSLS_soil_orgclay(1) && Z < ULSnSLS_soil_tergravel(1) 
                        Eb = ULSnSLS_soil_orgclay(8)*1000;
                        v = ULSnSLS_soil_orgclay(9);
                        G = Eb/(2*(1+v));
                    elseif Z >= ULSnSLS_soil_fill(1) && Z < ULSnSLS_soil_orgclay(1) 
                        Eb = ULSnSLS_soil_fill(8)*1000;   
                        v = ULSnSLS_soil_fill(9);
                        G = Eb/(2*(1+v));
                    end

                    %% Fa
                    mu = sqrt((2*pi*G)/(log(2*L/d)*Eb*A));
                    Kbi = (2*d*G)/(1-v);
                    omega = Kbi/(Eb*A*mu);
                    mu_L_2 = 2*mu*L; 
                    Fa = (mu_L_2 + sinh(mu_L_2)+(omega^2)*(sinh(mu_L_2)-mu_L_2)+2*omega*(cosh(mu_L_2)-1))/((2+2*omega^2)*sinh(mu_L_2)+4*omega*cosh(mu_L_2));
                    
                    %% aj
                    if isempty(find(pile_types==1)) == 1
                        inner_id = 0;
                    else
                        inner_id = find(pile_types==1);
                    end
                    if isempty(find(pile_types==2)) == 1
                        edge_id = 0;
                    else
                        edge_id = find(pile_types==2);
                    end
                    if isempty(find(pile_types==3)) == 1
                        corner_id = 0;
                    else
                        corner_id = find(pile_types==3);
                    end
                    list_of_piles_to_calc = [inner_id(1), edge_id(1), corner_id(1)];

                    aj_per_eq = [];
                    k = 1;
                    
                    for i_loops = 1:length(list_of_piles_to_calc)  % shows how many pile types are requried 
                        % skip analysis on same type of piles
                        if list_of_piles_to_calc(i_loops) == 0 && i_loops < 3
                            continue 
                        elseif list_of_piles_to_calc(i_loops) == 0 && i_loops == 3
                            break
                        end

                        % calculate aj for each pile types
                        aj_sum_each_pile_type = [0, 0, 0];
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
                    number_of_pile_types_to_analyse = (N_pile_I>0)+(N_pile_E>0)+(N_pile_C>0);
                    piles_to_design = unique(pile_types);
                    ULS_load_P_per_pile_type_P = [0; 0; 0];		% load distribution

                    if number_of_pile_types_to_analyse == 1
                        ULS_load_P_per_pile_type_P(piles_to_design) = ULS_load_total/N;
                        SLS_load_P_per_pile_type_P(piles_to_design) = SLS_load_total/N;
                        M = nonzeros(aj_per_eq);
                        Q_matrix = [SLS_load_total/N];

                    elseif number_of_pile_types_to_analyse == 2
                        M = [aj_per_eq(1:2,piles_to_design(1)),aj_per_eq(1:2,piles_to_design(2))];
                        ratio_of_loads = M\[1; 1];
                        if N_pile_I > 0 && N_pile_E > 0
                            n1 = N_pile_I;
                            nn1 = 1;
                            n2 = N_pile_E;
                            nn2 = 2;
                        elseif N_pile_I > 0 && N_pile_C > 0
                            n1 = N_pile_I;
                            nn1 = 1;
                            n2 = N_pile_C;
                            nn2 = 3;
                        elseif N_pile_E > 0 && N_pile_C > 0
                            n1 = N_pile_E;
                            nn1 = 2;
                            n2 = N_pile_C;
                            nn2 = 3;
                        end
                        Q1_ULS = ULS_load_total/([n1, n2]*ratio_of_loads);
                        ULS_load_P_ratio = Q1_ULS*ratio_of_loads;
                        ULS_load_P_per_pile_type_P(nn1) = ULS_load_P_ratio(1);
                        ULS_load_P_per_pile_type_P(nn2) = ULS_load_P_ratio(2);

                        Q1_SLS = SLS_load_total/([n1, n2]*ratio_of_loads);
                        SLS_load_P_ratio = Q1_SLS*ratio_of_loads;
                        SLS_load_P_per_pile_type_P(nn1) = SLS_load_P_ratio(1);
                        SLS_load_P_per_pile_type_P(nn2) = SLS_load_P_ratio(2);
                        Q_matrix = [SLS_load_P_ratio(1); SLS_load_P_ratio(2)];

                    elseif number_of_pile_types_to_analyse == 3
                        M = aj_per_eq;
                        ratio_of_loads = M\[1; 1; 1];
                        Q1_ULS = ULS_load_total/([N_pile_I, N_pile_E, N_pile_C]*ratio_of_loads);
                        ULS_load_P_per_pile_type_P = Q1_ULS*ratio_of_loads;

                        Q1_SLS = SLS_load_total/([N_pile_I, N_pile_E, N_pile_C]*ratio_of_loads);
                        Q_matrix = Q1_SLS*ratio_of_loads;
                        SLS_load_P_per_pile_type_P = Q_matrix;
                    end
                    
                    % calculations check
                    if round(ULS_load_total) ~= round(ULS_load_P_per_pile_type_P(1)*N_pile_I + ULS_load_P_per_pile_type_P(2)*N_pile_E + ULS_load_P_per_pile_type_P(3)*N_pile_C)
                        error('your sum of loads on each pile does not add up to total uls axial loads you are imposing')
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% loads on each pile types
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    ULS_load_pile_type_I = ULS_load_P_per_pile_type_P(1) + ULS_load_P_from_M_per_pile_type_I;
                    ULS_load_pile_type_E = ULS_load_P_per_pile_type_P(2) + ULS_load_P_from_M_per_pile_type_E;
                    ULS_load_pile_type_C = ULS_load_P_per_pile_type_P(3) + ULS_load_P_from_M_per_pile_type_C;
            		ULS_load_pile_type_M_only = [ULS_load_P_from_M_per_pile_type_I, ULS_load_P_from_M_per_pile_type_E, ULS_load_P_from_M_per_pile_type_C];
                    ULS_load_pile_type = [ULS_load_pile_type_I, ULS_load_pile_type_E, ULS_load_pile_type_C];

                    SLS_load_pile_type_I = SLS_load_P_per_pile_type_P(1) + SLS_load_P_from_M_per_pile_type_I;
                    SLS_load_pile_type_E = SLS_load_P_per_pile_type_P(2) + SLS_load_P_from_M_per_pile_type_E;
                    SLS_load_pile_type_C = SLS_load_P_per_pile_type_P(3) + SLS_load_P_from_M_per_pile_type_C;
                    SLS_load_pile_type_M_only = [SLS_load_P_from_M_per_pile_type_I, SLS_load_P_from_M_per_pile_type_E, SLS_load_P_from_M_per_pile_type_C];
                    SLS_load_pile_type = [SLS_load_pile_type_I, SLS_load_pile_type_E, SLS_load_pile_type_C];

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% base resistance - group piles
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    % base on terrace gravel
                    qb_group = 0;
                    if ULSnSLS_soil_Lonclay1(1) > Z && Z >= ULSnSLS_soil_tergravel(1)
                        Nq = exp(pi*tand(ULSnSLS_soil_tergravel(5)))*(tand(45+0.5*ULSnSLS_soil_tergravel(5))).^2;
                        Ngam = 2*(Nq-1)*tand(ULSnSLS_soil_tergravel(5));
                        % overburden stress
                        overburden_stress = (d-ULSnSLS_soil_tergravel(1))*ULSnSLS_soil_tergravel(3)+ ...
                            ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+ (ULSnSLS_soil_fill(1)-basement_depth)*ULSnSLS_soil_fill(3);
                        qb_group = Nq*overburden_stress+0.4*Ngam*ULSnSLS_soil_tergravel(3)*bb;
                    % base on London clay
                    elseif ULSnSLS_soil_Lonclay2(1) > Z && Z >= ULSnSLS_soil_Lonclay1(1)
                        if (1+0.2*(L/bb)) > 1.5
                            square_part = 1.5;
                        else
                            square_part = (1+0.2*(L/bb));
                        end
                        Nc = square_part*5*(1+0.2*(ll/bb));
                        qb_group = Nc*ULSnSLS_soil_Lonclay1(4);
                    % base on stiff London clay
                    elseif ULSnSLS_soil_Lonclay3(1) > Z && Z >= ULSnSLS_soil_Lonclay2(1)
                        if (1+0.2*(L/bb)) > 1.5
                            square_part = 1.5;
                        else
                            square_part = (1+0.2*(L/bb));
                        end
                        Nc = square_part*5*(1+0.2*(ll/bb));
                        qb_group = Nc*ULSnSLS_soil_Lonclay2(4);
                    % base on claystone
                    elseif (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2)) > Z && Z >= ULSnSLS_soil_Lonclay3(1)
                        if (1+0.2*(L/bb)) > 1.5
                            square_part = 1.5;
                        else
                            square_part = (1+0.2*(L/bb));
                        end
                        Nc = square_part*5*(1+0.2*(ll/bb));
                        qb_group = Nc*ULSnSLS_soil_Lonclay3(4);
                    end

                    % static base resistance force
                    ULS_Rb_group = Ab*qb_group;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% shaft resistance - group pile
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    % shaft resistance depth
                    t_shaft_dep = Z;

                    %% shaft on fill
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_fill(1) && t_shaft_dep < ULSnSLS_soil_orgclay(1)
                        len_fill = t_shaft_dep - basement_depth;
                    elseif t_shaft_dep >= ULSnSLS_soil_fill(1) && t_shaft_dep >= ULSnSLS_soil_orgclay(1)
                        len_fill = ULSnSLS_soil_fill(1) + ULSnSLS_soil_fill(2) - basement_depth;
                    else
                        len_fill = 0;
                    end
                    % beta
                    beta_fill = Ks*tand(del_fill);
                    % pressure
                    qs_fill = beta_fill*(len_fill*0.5)*ULSnSLS_soil_fill(3);
                    % force
                    Qs_fill = qs_fill*perimeter*len_fill;

                    %% shaft on organic clay
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_orgclay(1) && t_shaft_dep < ULSnSLS_soil_tergravel(1)
                        len_orgclay = t_shaft_dep-ULSnSLS_soil_orgclay(1);
                    elseif t_shaft_dep >= ULSnSLS_soil_orgclay(1) && t_shaft_dep >= ULSnSLS_soil_tergravel(1)
                        len_orgclay = ULSnSLS_soil_orgclay(2);
                    else
                        len_orgclay = 0;
                    end    
                    % pressure
                    if ULSnSLS_soil_orgclay(7)*ULSnSLS_soil_orgclay(4) > 110
                        qs_org_clay = 110;
                    else
                        qs_org_clay = ULSnSLS_soil_orgclay(7)*ULSnSLS_soil_orgclay(4);
                    end
                    % Force
                    Qs_org_clay = qs_org_clay*perimeter*len_orgclay;

                    %% shaft on terrace gravel
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_tergravel(1) && t_shaft_dep < ULSnSLS_soil_Lonclay1(1)
                        len_tergravel = t_shaft_dep-ULSnSLS_soil_tergravel(1);
                    elseif t_shaft_dep >= ULSnSLS_soil_tergravel(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay1(1)
                        len_tergravel = ULSnSLS_soil_tergravel(2);
                    else
                        len_tergravel = 0;
                    end    
                    % beta
                    beta_tergravel = Ks*tand(del_tergravel);
                    % pressure
                    overburden_pres = (len_tergravel*0.5)*ULSnSLS_soil_tergravel(3)+...
                                       ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+...
                                      (ULSnSLS_soil_fill(1)+ULSnSLS_soil_fill(2)-basement_depth)*ULSnSLS_soil_fill(3);
                    qs_tergravel = beta_tergravel*overburden_pres;
                    % force
                    Qs_tergravel = qs_tergravel*perimeter*len_tergravel;

                    %% shaft on London clay
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep < ULSnSLS_soil_Lonclay2(1)
                        len_Lonclay1 = t_shaft_dep-ULSnSLS_soil_Lonclay1(1);
                    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay2(1)
                        len_Lonclay1 = ULSnSLS_soil_Lonclay1(2);
                    else
                        len_Lonclay1 = 0;
                    end    
                    % pressure
                    if ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4) > 110
                        qs_Lonclay1 = 110;
                    else
                        qs_Lonclay1 = ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4);
                    end
                    % force
                    Qs_Lonclay1 = qs_Lonclay1*perimeter*len_Lonclay1;

                    %% shaft on Stiff London clay
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep < ULSnSLS_soil_Lonclay3(1)
                        len_Lonclay2 = t_shaft_dep-ULSnSLS_soil_Lonclay2(1);
                    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay3(1)
                        len_Lonclay2 = ULSnSLS_soil_Lonclay2(2);
                    else
                        len_Lonclay2 = 0;
                    end    
                    % pressure
                    if ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4) > 110
                        qs_Lonclay2 = 110;
                    else
                        qs_Lonclay2 = ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4);
                    end
                    % force
                    Qs_Lonclay2 = qs_Lonclay2*perimeter*len_Lonclay2;

                    %% shaft on Claystone
                    % length
                    if t_shaft_dep >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep < (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
                        len_Lonclay3 = t_shaft_dep-ULSnSLS_soil_Lonclay3(1);
                    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep >= (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
                        len_Lonclay3 = ULSnSLS_soil_Lonclay3(2);
                    else
                        len_Lonclay3 = 0;
                    end 
                    % pressure
                    if ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4) > 110
                        qs_Lonclay3 = 110;
                    else
                        qs_Lonclay3 = ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4);
                    end
                    % force
                    Qs_Lonclay3 = qs_Lonclay3*perimeter*len_Lonclay3;

                    %% static base resistance force
                    ULS_Rs_group = -Qs_fill-Qs_org_clay+Qs_tergravel+Qs_Lonclay1+Qs_Lonclay2+Qs_Lonclay3;
                    % fill and organic clay layer causes negative skin friction

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Ultimate resistance - group piles
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    % compression
                    % ULS_Ru_from_Rs_Rb = ULS_Rb_group/R_EQU_base + ULS_Rs_group/R_EQU_sh_comp;    % unit: kN
                    % ULS_Ru_from_Rs = ULS_Rs_group/1.2;                               % unit: kN
                    % ULS_Ru_from_fck_cube = A*fck_cube*0.25;                    % unit: kN
                    % ULS_Ru_comp = min([ULS_Ru_from_Rs_Rb;ULS_Ru_from_Rs;ULS_Ru_from_fck_cube]);
                    ULS_Ru_comp_group = ULS_Rb_group/STR_GEO_AMR{comb{i}(3)}(1) + ULS_Rs_group/STR_GEO_AMR{comb{i}(3)}(2);    % unit: kN
                    % tension
                    ULS_Ru_ten_group = ULS_Rs_group/STR_GEO_AMR{comb{i}(3)}(3);
                    ULS_Ru_group = [ULS_Ru_comp_group, ULS_Ru_ten_group];

        			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		    %% base resistance stress - single
        		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		    qb = 0;
        		    if ULSnSLS_soil_Lonclay1(1) > Z && Z >= ULSnSLS_soil_tergravel(1)
        		        Nq = exp(pi*tand(ULSnSLS_soil_tergravel(5)))*(tand(45+0.5*ULSnSLS_soil_tergravel(5))).^2;
                        Nq_SLS = exp(pi*tand(soil(3,5)))*(tand(45+0.5*soil(3,5))).^2;
        		        % overburden stress
        		        overburden_stress = (Z-ULSnSLS_soil_tergravel(1))*ULSnSLS_soil_tergravel(3)+ ...
        		            ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+ (ULSnSLS_soil_fill(1)-basement_depth)*ULSnSLS_soil_fill(3);
                        overburden_stress_SLS = (Z-soil(3,1))*soil(3,3)+soil(2,2)*soil(2,3)+ (soil(1,1)-basement_depth)*soil(1,3);
        		        qb = Nq*overburden_stress;
                        qb_SLS = Nq_SLS*overburden_stress_SLS;
        		    % base on London clay
        		    elseif ULSnSLS_soil_Lonclay2(1) > Z && Z >= ULSnSLS_soil_Lonclay1(1)
        		        qb = 9*ULSnSLS_soil_Lonclay1(4);
                        qb_SLS = 9*soil(4,4);
        		    % base on stiff London clay
        		    elseif ULSnSLS_soil_Lonclay3(1) > Z && Z >= ULSnSLS_soil_Lonclay2(1)
        		        qb = 9*ULSnSLS_soil_Lonclay2(4);
                        qb_SLS = 9*soil(5,4);
        		    % base on claystone
        		    elseif (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2)) > Z && Z >= ULSnSLS_soil_Lonclay3(1)
        		        qb = 9*ULSnSLS_soil_Lonclay3(4);
                        qb_SLS = 9*soil(5,4);
        		    end
        		    % static base resistance force
        		    ULS_Rb_single = A*qb;
                    SLS_Rb_single = A*qb_SLS;

        		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		    %% shaft resistance stress
        		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		    t_shaft_dep = L;	    % shaft resistance depth
        		    
        		    %% shaft on fill
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_fill(1) && t_shaft_dep < ULSnSLS_soil_orgclay(1)
        		        len_fill = t_shaft_dep - basement_depth;
        		    elseif t_shaft_dep >= ULSnSLS_soil_fill(1) && t_shaft_dep >= ULSnSLS_soil_orgclay(1)
        		        len_fill = ULSnSLS_soil_fill(1)+ULSnSLS_soil_fill(2) - basement_depth;
        		    else
        		        len_fill = 0;
        		    end
        		    % beta
        		    beta_fill = Ks*tand(del_fill);
        		    % pressure
        		    qs_fill = beta_fill*(len_fill*0.5)*ULSnSLS_soil_fill(3);
                    qs_fill_SLS = beta_fill*(len_fill*0.5)*soil(1,3);
        		    % force
        		    Qs_fill = qs_fill*pi*d*len_fill;
                    Qs_fill_SLS = qs_fill_SLS*pi*d*len_fill;

        		    %% shaft on organic clay
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_orgclay(1) && t_shaft_dep < ULSnSLS_soil_tergravel(1)
        		        len_orgclay = t_shaft_dep-ULSnSLS_soil_orgclay(1);
        		    elseif t_shaft_dep >= ULSnSLS_soil_orgclay(1) && t_shaft_dep >= ULSnSLS_soil_tergravel(1)
        		        len_orgclay = ULSnSLS_soil_orgclay(2);
        		    else
        		        len_orgclay = 0;
        		    end    
        		    % pressure
        		    if ULSnSLS_soil_orgclay(7)*ULSnSLS_soil_orgclay(4) > 110
        		        qs_org_clay = 110;
        		    else
        		        qs_org_clay = ULSnSLS_soil_orgclay(7)*ULSnSLS_soil_orgclay(4);
        		    end
                    if soil(2,7)*soil(2,4) > 110
                        qs_org_clay_SLS = 110;
                    else
                        qs_org_clay_SLS = soil(2,7)*soil(2,4);
                    end
        		    % force
        		    Qs_org_clay = qs_org_clay*pi*d*len_orgclay;
                    Qs_org_clay_SLS = qs_org_clay_SLS*pi*d*len_orgclay;

        		    %% shaft on terrace gravel
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_tergravel(1) && t_shaft_dep < ULSnSLS_soil_Lonclay1(1)
        		        len_tergravel = t_shaft_dep-ULSnSLS_soil_tergravel(1);
        		    elseif t_shaft_dep >= ULSnSLS_soil_tergravel(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay1(1)
        		        len_tergravel = ULSnSLS_soil_tergravel(2);
        		    else
        		        len_tergravel = 0;
        		    end    
        		    % beta
        		    beta_tergravel = Ks*tand(del_tergravel);
        		    % pressure
        		    overburden_pres = (len_tergravel*0.5)*ULSnSLS_soil_tergravel(3)+...
        		                       ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+...
        		                      (ULSnSLS_soil_fill(1)+ULSnSLS_soil_fill(2)-basement_depth)*ULSnSLS_soil_fill(3);
        		    qs_tergravel = beta_tergravel*overburden_pres;
                    overburden_pres_SLS = (len_tergravel*0.5)*soil(3,3)+soil(2,2)*soil(2,3)+(soil(1,1)+soil(1,2)-basement_depth)*soil(1,3);
                    qs_tergravel_SLS = beta_tergravel*overburden_pres_SLS;
        		    % force
        		    Qs_tergravel = qs_tergravel*pi*d*len_tergravel;
        		    Qs_tergravel_SLS = qs_tergravel_SLS*pi*d*len_tergravel;

        		    %% shaft on London clay
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep < ULSnSLS_soil_Lonclay2(1)
        		        len_Lonclay1 = t_shaft_dep-ULSnSLS_soil_Lonclay1(1);
        		    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay2(1)
        		        len_Lonclay1 = ULSnSLS_soil_Lonclay1(2);
        		    else
        		        len_Lonclay1 = 0;
        		    end    
        		    % pressure
        		    if ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4) > 110
        		        qs_Lonclay1 = 110;
        		    else
        		        qs_Lonclay1= ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4);
        		    end
                    if soil(4,7)*soil(4,4) > 110
                        qs_Lonclay1_SLS = 110;
                    else
                        qs_Lonclay1_SLS = soil(4,7)*soil(4,4) ;
                    end
        		    % force
        		    Qs_Lonclay1 = qs_Lonclay1*pi*d*len_Lonclay1;
                    Qs_Lonclay1_SLS = qs_Lonclay1_SLS*pi*d*len_Lonclay1;

        		    %% shaft on Stiff London clay
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep < ULSnSLS_soil_Lonclay3(1)
        		        len_Lonclay2 = t_shaft_dep-ULSnSLS_soil_Lonclay2(1);
        		    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep >= ULSnSLS_soil_Lonclay3(1)
        		        len_Lonclay2 = ULSnSLS_soil_Lonclay2(2);
        		    else
        		        len_Lonclay2 = 0;
        		    end    
        		    % pressure
        		    if ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4) > 110
        		        qs_Lonclay2 = 110;
        		    else
        		        qs_Lonclay2= ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4);
        		    end
                    if soil(5,7)*soil(5,4) > 110
                        qs_Lonclay2_SLS = 110;
                    else
                        qs_Lonclay2_SLS = soil(5,7)*soil(5,4);
                    end
        		    % force
        		    Qs_Lonclay2 = qs_Lonclay2*pi*d*len_Lonclay2;
                    Qs_Lonclay2_SLS = qs_Lonclay2_SLS*pi*d*len_Lonclay2;

        		    %% shaft on Claystone
        		    % length
        		    if t_shaft_dep >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep < (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
        		        len_Lonclay3 = t_shaft_dep-ULSnSLS_soil_Lonclay3(1);
        		    elseif t_shaft_dep >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep >= (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
        		        len_Lonclay3 = ULSnSLS_soil_Lonclay3(2);
        		    else
        		        len_Lonclay3 = 0;
        		    end 
        		    % pressure
        		    if ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4) > 110
        		        qs_Lonclay3 = 110;
        		    else
        		        qs_Lonclay3= ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4);
        		    end
                    if soil(6,7)*soil(6,4) > 110
                        qs_Lonclay3_SLS = 110;
                    else
                        qs_Lonclay3_SLS = soil(6,7)*soil(6,4);
                    end
        		    % force
        		    Qs_Lonclay3 = qs_Lonclay3*pi*d*len_Lonclay3;
                    Qs_Lonclay3_SLS = qs_Lonclay3_SLS*pi*d*len_Lonclay3;
        		    
        		    %% static shaft resistance force
        		    ULS_Rs_single = -Qs_fill-Qs_org_clay+Qs_tergravel+Qs_Lonclay1+Qs_Lonclay2+Qs_Lonclay3;
                    SLS_Rs_single = -Qs_fill_SLS-Qs_org_clay_SLS+Qs_tergravel_SLS+Qs_Lonclay1_SLS+Qs_Lonclay2_SLS+Qs_Lonclay3_SLS;
        		    % fill and organic clay layer causes negative skin friction

        		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		    %% Ultimate resistance - individual piles
        		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        			% compression bearing capacity	(LSDA Foundation No 1, 3.2)
        		    ULS_single_from_Rs_Rb = ULS_Rb_single/STR_GEO_AMR{comb{i}(3)}(1) + ULS_Rs_single/STR_GEO_AMR{comb{i}(3)}(2);    % unit: kN
        		    ULS_single_from_Rs = ULS_Rs_single/1.2;                               % unit: kN
        		    ULS_single_from_fck_cube = A*fck_cube*0.25;                    % unit: kN
        		    ULS_Ru_comp_single = min([ULS_single_from_Rs_Rb;ULS_single_from_Rs;ULS_single_from_fck_cube]);
        	    	ULS_Ru_ten_single = ULS_Rs_single/STR_GEO_AMR{comb{i}(3)}(3);

                    ULS_Ru_single = [ULS_Ru_comp_single, ULS_Ru_ten_single];
                    SLS_Ru_single = [SLS_Rb_single+SLS_Rs_single, SLS_Rs_single];

        	    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	        %% Settlement - individual piles
        	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	        pile_types_settlement_from_moments = [0, 0, 0];				% unit: mm
        	        Kpile = [0, 0, 0];											% unit: Kn/m
        	    	for pile_types_settlement = 1:2*length(pile_types_settlement_from_moments)

        	    		% load on each pile type
        	    		if pile_types_settlement <= 3
        	    			Q = SLS_load_pile_type_M_only(pile_types_settlement);
        	    		elseif pile_types_settlement > 3 && pile_types_settlement <= 2*length(pile_types_settlement_from_moments)
        	    			if SLS_load_total > 0
                                Q = 100;
                            elseif SLS_load_total < 0  
                                Q = -100;
                            end
                        end
                        
                        if Q == 0 && pile_types_settlement < 2*length(pile_types_settlement_from_moments)
                            continue
                        elseif Q == 0 && pile_types_settlement < 2*length(pile_types_settlement_from_moments)
                            break
                        end
                        
                        if SLS_load_total > 0
                            SLS_Ru = SLS_Ru_single(1);
                        elseif SLS_load_total < 0  
                            SLS_Ru = SLS_Ru_single(2);
                        end
                        
        		    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		        %% Settlement - individual piles
        		        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		        % soil under base Elasticity
        			    Eb = 0;
        			    if Z>= ULSnSLS_soil_Lonclay3(1) && Z< ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2) 
        			        Eb = ULSnSLS_soil_Lonclay3(8)*1000;
        			    elseif Z>= ULSnSLS_soil_Lonclay2(1) && Z< ULSnSLS_soil_Lonclay3(1) 
        			        Eb = ULSnSLS_soil_Lonclay2(8)*1000;
        			    elseif Z>= ULSnSLS_soil_Lonclay1(1) && Z< ULSnSLS_soil_Lonclay2(1) 
        			        Eb = ULSnSLS_soil_Lonclay1(8)*1000;
        			    elseif Z>= ULSnSLS_soil_tergravel(1) && Z< ULSnSLS_soil_Lonclay1(1) 
        			        Eb = ULSnSLS_soil_tergravel(8)*1000;
        			    elseif Z>= ULSnSLS_soil_orgclay(1) && Z< ULSnSLS_soil_tergravel(1) 
        			        Eb = ULSnSLS_soil_orgclay(8)*1000;
        			    elseif Z>= ULSnSLS_soil_fill(1) && Z< ULSnSLS_soil_orgclay(1) 
        			        Eb = ULSnSLS_soil_fill(8)*1000;     
        			    end

        		        % parameters (5 greek)
        		        alpha = SLS_Rs_single;
        		        if Q > 0
        		        	beta = SLS_Rb_single*d*Eb;
        		        	delta = SLS_Rb_single*0.6;
        		        elseif Q < 0
        		        	beta = 0;
        		        	delta = 0;
        		        end
        		        lambda = d*Ms;
        		        eta = d*Eb;
        		        
        		        % parameters (3 alphAet)
        		        aa = eta*(Q-alpha)-beta;
        		        bb = Q*(delta+lambda*eta)-alpha*delta-beta*lambda;
        		        cc = lambda*delta*Q;
        		        
        		        % soil settlement (sr)
        		        if abs(Q) <= abs(SLS_Ru)
        		            srp = 1000*((-1*bb)+sqrt((bb.^2)-(4*aa*cc)))/(2*aa);
        		            sre = 1000*((-1*bb)-sqrt((bb.^2)-(4*aa*cc)))/(2*aa);
        		            if Q > 0
                                if sre > 0 && srp <= 0
                                    sr = sre;
                                elseif srp > 0 && sre <= 0
                                    sr = srp;
                                else
                                    sr = 9999;
                                end
                            elseif Q < 0
                                if sre < 0 && srp >= 0
                                    sr = sre;
                                elseif srp < 0 && sre >= 0
                                    sr = srp;
                                else
                                    sr = 9999;
                                end
                            end
        		        elseif abs(Q) > abs(SLS_Ru)
        		            sr = 9999;
        		        end
        		        
        		        % pile compression elastic settlement (se)
        		        if abs(Q) <= abs(SLS_Rs_single)
        		            se = 1000*(4*Q*(Lo + Ke*(L-Lo)))/(pi*Ep*(d.^2));
        		        elseif abs(Q) > abs(ULS_Rs_single) && abs(Q) <= abs(SLS_Ru)
        		            se = 1000*(4*((Q*L)-((ULS_Rs_single)*(1-Ke)*(L-Lo))))/(pi*Ep*(d.^2));
        		        elseif abs(Q) > abs(SLS_Ru)
        		            se = 9999;
        		        end
        		        
        		        % modify datalog
        		        if se==9999 || sr==9999
        		            s_single = 9999;
        		        else
        		            s_single = sr + se;
        		        end

        		        % collect results
        		        % from moment
        		        if pile_types_settlement <= 3
        		        	if Q == 0
        		        		pile_types_settlement_from_moments(pile_types_settlement) = 0;
        		        	else
        	    				pile_types_settlement_from_moments(pile_types_settlement) = s_single;
                            end
        	    		% to find the Kpile
        	    		elseif pile_types_settlement > 3 && pile_types_settlement <= 2*length(pile_types_settlement_from_moments)
        	    			Kpile(pile_types_settlement-3) = Q/(s_single/1000);
                        end
        		    end

        		    % calculation check 
        		    if not(isequal(Kpile(1),Kpile(2),Kpile(3)))
        		    	error('Kpile should be equal for each piles')
        		    elseif isequal(Kpile(1),Kpile(2),Kpile(3))
        		    	Kpile = Kpile(1);
        		    end

        	        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% settlement - group pile effect - based on Randolpy 1994 / From Craig pp 354-355
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    settlements_from_P = (1/Kpile)*M*Q_matrix;			% unit: m
                    pile_types_settlements_from_P_group = [0, 0, 0];
                    % calculation check and arrange to each pile types
                    if length(settlements_from_P) == 1
                        if N_pile_I > 0 
                            nn1 = 1;
                        elseif N_pile_E > 0
                            nn1 = 2;
                        elseif N_pile_C > 0
                            nn1 = 3;
                        end
                    	pile_types_settlements_from_P_group(nn1) = settlements_from_P(1);
                    elseif length(settlements_from_P) == 2
        			    if round(settlements_from_P(1),3,'significant') ~= round(settlements_from_P(2),3,'significant')
        			    	error('Kpile should be equal for each piles')
        			    elseif round(settlements_from_P(1),3,'significant') == round(settlements_from_P(2),3,'significant')
        			    	if N_pile_I > 0 && N_pile_E > 0
                                nn1 = 1;
                                nn2 = 2;
                            elseif N_pile_I > 0 && N_pile_C > 0
                                nn1 = 1;
                                nn2 = 3;
                            elseif N_pile_E > 0 && N_pile_C > 0
                                nn1 = 2;
                                nn2 = 3;
                            end
                            pile_types_settlements_from_P_group(nn1) = settlements_from_P(1);
                            pile_types_settlements_from_P_group(nn2) = settlements_from_P(2);
        			    end
        			elseif length(settlements_from_P) == 3
        			    if not(isequal(round(settlements_from_P(1),3,'significant'),round(settlements_from_P(2),3,'significant'),round(settlements_from_P(3),3,'significant')))
        			    	disp(settlements_from_P)
                            error('Kpile should be equal for each piles')
        			    elseif isequal(round(settlements_from_P(1),3,'significant'),round(settlements_from_P(2),3,'significant'),round(settlements_from_P(3),3,'significant'))
        			    	pile_types_settlements_from_P_group = settlements_from_P';
        			    end
        			end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% overall settlement - group pile effect and from additional load from moments
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    overall_settlement_for_each_pile = pile_types_settlement_from_moments + 1000*abs(pile_types_settlements_from_P_group);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% bearing & settlement check and tabulations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % slim_group - a number that limit maximum settlement
                    %(slim_percent_group/100)*d*1000 - a percentage of the pile diameters


                    % ULS group 
                    if abs(ULS_Ru_group(loads_cases)) >= abs(ULS_load_total) 
                        % ULS and SLS single
                        if N_pile_I == 0 
                            inner_pile_ok = 9999;
                        elseif N_pile_I > 0 
                            if abs(ULS_Ru_single(loads_cases)) >= ULS_load_pile_type(1) 
                                if (slim_percent_group/100)*d*1000 < overall_settlement_for_each_pile(1) || overall_settlement_for_each_pile(1) == 0  % unit: mm 
                    	            inner_pile_ok = 0;
                                elseif (slim_percent_group/100)*d*1000 >= overall_settlement_for_each_pile(1)	% unit: mm 
                                    inner_pile_ok = 1;
                                end
                            else 
                                inner_pile_ok = 0;
                            end
                        end
                        
                        if N_pile_E == 0
                            edge_pile_ok = 9999;
                        elseif N_pile_E > 0
                            if abs(ULS_Ru_single(loads_cases)) >= ULS_load_pile_type(2) 
                    	        if (slim_percent_group/100)*d*1000 < overall_settlement_for_each_pile(2) || overall_settlement_for_each_pile(2) == 0  % unit: mm 
                    	            edge_pile_ok = 0;
                                elseif (slim_percent_group/100)*d*1000 >= overall_settlement_for_each_pile(2) 
                                    edge_pile_ok = 1;
                                end
                            else 
                                edge_pile_ok = 0;
                            end
                        end
                        
                        if N_pile_C == 0
                            corner_pile_ok = 9999;
                        elseif N_pile_C > 0
                            if abs(ULS_Ru_single(loads_cases)) >= ULS_load_pile_type(3)  
                    	        if (slim_percent_group/100)*d*1000 < overall_settlement_for_each_pile(3) || overall_settlement_for_each_pile(3) == 0  % unit: mm 
                            	    corner_pile_ok = 0;
                                elseif (slim_percent_group/100)*d*1000 >= overall_settlement_for_each_pile(3) 
                                    corner_pile_ok = 1;
                                end
                            else 
                                corner_pile_ok = 0;
                            end
                        end           
                    else 
                        inner_pile_ok = 0;
                        edge_pile_ok = 0;
                        corner_pile_ok = 0;
                    end
                    % tabulated every calculations results for numerical checks
                    dimensions = [Z, L, d, N_pile_I, N_pile_E, N_pile_C, N, A, dx, bb, ll, Ab, perimeter, As, Ixx, Iyy, orientation]; 
                    loads = [critical_load_combi, ULS_load_total, ULS_load_pile_type_M_only, ULS_load_P_per_pile_type_P(1), ULS_load_P_per_pile_type_P(2), ULS_load_P_per_pile_type_P(3), ...
                            ULS_load_pile_type, SLS_load_total, SLS_load_pile_type_M_only, SLS_load_P_per_pile_type_P(1), SLS_load_P_per_pile_type_P(2), SLS_load_P_per_pile_type_P(3), SLS_load_pile_type];
                    resistance = [ULS_Rb_group, ULS_Rs_group, ULS_Ru_group(1), -ULS_Ru_group(2), ULS_Rb_single, ULS_Rs_single, ULS_single_from_Rs_Rb, ULS_single_from_Rs, ULS_single_from_fck_cube, ULS_Ru_single(1), -ULS_Ru_single(2)];
                    settlement = [pile_types_settlement_from_moments, pile_types_settlements_from_P_group, overall_settlement_for_each_pile]; 
                    pile_check = [inner_pile_ok, edge_pile_ok, corner_pile_ok];
                    cost_of_foundation = [L*0.25*pi*d.^2, N*L*0.25*pi*d.^2];

                    every_results(tab,:) = [dimensions, loads, resistance, settlement, pile_check, cost_of_foundation];
                    tab = tab+1;

                    if inner_pile_ok == 1 || edge_pile_ok == 1 || corner_pile_ok == 1
        	            if i == 1  % combination 1
        		            if ULS_load_total > 0
        		                ULSnSLS_design_C1_comp(row,:) = [dimensions, loads, resistance, settlement, pile_check, cost_of_foundation];
        		                row = row+1;
        	                elseif ULS_load_total < 0
        		                ULSnSLS_design_C1_ten(row,:) = [dimensions, loads, resistance, settlement, pile_check, cost_of_foundation];
        		                row = row+1;
        		            end
        	            elseif i == 2  % combination 2
        		            if ULS_load_total > 0
        		                ULSnSLS_design_C2_comp(row,:) = [dimensions, loads, resistance, settlement, pile_check, cost_of_foundation];
        		                row = row+1;
        	                elseif ULS_load_total < 0
        		                ULSnSLS_design_C2_ten(row,:) = [dimensions, loads, resistance, settlement, pile_check, cost_of_foundation];
        		                row = row+1;
        	                end
                        end
        	        else
        	            continue
        	        end
                end
            end
	    end
    end
end
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the worst combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loads_cases = 1:2
    kk = 1;
    n = 1;
    %% compare settlement results
    if ULS_Load_total_list(loads_cases) > 0 && (isempty(ULSnSLS_design_C1_comp) == 0 && isempty(ULSnSLS_design_C2_comp) == 0)
        [w1,~]=size(ULSnSLS_design_C1_comp);
        [w2,~]=size(ULSnSLS_design_C2_comp);
        % these variables stores the number of times a combination gives a worse scenario
        for l = 1:w1
            for m = n:w2
                if ULSnSLS_design_C1_comp(l,3) == ULSnSLS_design_C2_comp(m,3) && ULSnSLS_design_C1_comp(l,7) == ULSnSLS_design_C2_comp(m,7) && ULSnSLS_design_C1_comp(l,2) == ULSnSLS_design_C2_comp(m,2)
                    if ULSnSLS_design_C1_comp(l,41) >= ULSnSLS_design_C2_comp(m,41) && ULSnSLS_design_C1_comp(l,48) >= ULSnSLS_design_C2_comp(m,48) 
                        ULSnSLS_Design_comp(kk,:) = [ULSnSLS_design_C2_comp(m,:), 2]; 
                        kk = kk + 1;
                        n = m;
                    elseif ULSnSLS_design_C1_comp(l,41) < ULSnSLS_design_C2_comp(m,41) && ULSnSLS_design_C1_comp(l,48) < ULSnSLS_design_C2_comp(m,48) 
                        ULSnSLS_Design_comp(kk,:) = [ULSnSLS_design_C1_comp(l,:), 1];
                        kk = kk + 1;
                        n = m;
                    end
                end
            end
        end 

    elseif ULS_Load_total_list(loads_cases) < 0 && (isempty(ULSnSLS_design_C1_ten) == 0 && isempty(ULSnSLS_design_C2_ten) == 0)
        [w1,~]=size(ULSnSLS_design_C1_ten);
        [w2,~]=size(ULSnSLS_design_C2_ten);
        % these variables stores the number of times a combination gives a worse scenario
        for l = 1:w1
            for m = n:w2
                if ULSnSLS_design_C1_ten(l,3) == ULSnSLS_design_C2_ten(m,3) && ULSnSLS_design_C1_ten(l,7) == ULSnSLS_design_C2_ten(m,7) && ULSnSLS_design_C1_ten(l,2) == ULSnSLS_design_C2_ten(m,2)
                    if abs(ULSnSLS_design_C1_ten(l,42)) >= abs(ULSnSLS_design_C2_ten(m,42)) &&  abs(ULSnSLS_design_C1_ten(l,49)) >= abs(ULSnSLS_design_C2_ten(m,49)) 
                        ULSnSLS_Design_ten(kk,:) = [ULSnSLS_design_C2_ten(m,:), 2]; 
                        kk = kk + 1;
                        n = m;
                    elseif abs(ULSnSLS_design_C1_ten(l,42)) < abs(ULSnSLS_design_C2_ten(m,42)) &&  abs(ULSnSLS_design_C1_ten(l,49)) < abs(ULSnSLS_design_C2_ten(m,49)) 
                        ULSnSLS_Design_ten(kk,:) = [ULSnSLS_design_C1_ten(l,:), 1];
                        kk = kk + 1;
                        n = m;
                    end
                end
            end
        end 

    else
        if loads_cases == 1
            ULSnSLS_design_C1_comp = [9999];
            ULSnSLS_design_C2_comp = [9999];
            ULSnSLS_Design_comp = [9999];
            %sorted_ULSnSLS_Design_comp = [9999];
            %ULSnSLS_Group_Design_comp = [9999];
            %sorted_ULSnSLS_Group_Design_comp = [9999];
            continue
        elseif loads_cases == 2
            ULSnSLS_design_C1_ten = [9999];
            ULSnSLS_design_C2_ten = [9999];
            ULSnSLS_Design_ten = [9999];
            %sorted_ULSnSLS_Design_ten = [9999];
            %ULSnSLS_Group_Design_ten = [9999];
            %sorted_ULSnSLS_Group_Design_ten = [9999];
            break
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sorted by number of piles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ULSnSLS_Design_comp(1,1) ~= 9999 
    %% number of piles
    pile_number_ULSnSLS_Design_comp = ULSnSLS_Design_comp;
    [~,order] = sort(pile_number_ULSnSLS_Design_comp(:,7));        
    pile_number_ULSnSLS_Design_comp = pile_number_ULSnSLS_Design_comp(order,:);
else 
    pile_number_ULSnSLS_Design_comp = [9999];
end

if ULSnSLS_Design_ten(1,1) ~= 9999 
    %% number of piles
    pile_number_ULSnSLS_Design_ten = ULSnSLS_Design_ten;
    [~,order] = sort(pile_number_ULSnSLS_Design_ten(:,7));        
    pile_number_ULSnSLS_Design_ten = pile_number_ULSnSLS_Design_ten(order,:);
else 
    pile_number_ULSnSLS_Design_ten = [9999];
end

%% cost analysis according to concrete volume used
if ULSnSLS_Design_comp(1,1) ~= 9999 
    [~,cost_col] = size(ULSnSLS_Design_comp);
    % individual piles - compression
    per_pile_cost_ULSnSLS_Design_comp = ULSnSLS_Design_comp;
    [~,order] = sort(per_pile_cost_ULSnSLS_Design_comp(:,cost_col-2));        
    per_pile_cost_ULSnSLS_Design_comp = per_pile_cost_ULSnSLS_Design_comp(order,:);
    % group piles - compression
    group_pile_cost_ULSnSLS_Design_comp = ULSnSLS_Design_comp;
    [~,order] = sort(group_pile_cost_ULSnSLS_Design_comp(:,cost_col-1));        
    group_pile_cost_ULSnSLS_Design_comp = group_pile_cost_ULSnSLS_Design_comp(order,:);
else 
    per_pile_cost_ULSnSLS_Design_comp = [9999];
    group_pile_cost_ULSnSLS_Design_comp = [9999];
end

if ULSnSLS_Design_ten(1,1) ~= 9999
    [~,cost_col] = size(ULSnSLS_Design_ten);
    % individual piles - tension
    per_pile_cost_ULSnSLS_Design_ten = ULSnSLS_Design_ten;
    [~,order] = sort(per_pile_cost_ULSnSLS_Design_ten(:,cost_col-2));        
    per_pile_cost_ULSnSLS_Design_ten = per_pile_cost_ULSnSLS_Design_ten(order,:);
    % group piles - tension
    group_pile_cost_ULSnSLS_Design_ten = ULSnSLS_Design_ten;
    [~,order] = sort(group_pile_cost_ULSnSLS_Design_ten(:,cost_col-1));        
    group_pile_cost_ULSnSLS_Design_ten = group_pile_cost_ULSnSLS_Design_ten(order,:);
else 
    per_pile_cost_ULSnSLS_Design_ten = [9999];
    group_pile_cost_ULSnSLS_Design_ten = [9999];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ULSnSLS calculation summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ULSnSLS_Group_Design = {ULSnSLS_design_C1_comp, ULSnSLS_design_C2_comp, ULSnSLS_design_C1_ten, ULSnSLS_design_C2_ten, ULSnSLS_Design_comp, ULSnSLS_Design_ten};
%ULSnSLS_Group_Design = {ULSnSLS_Group_Design_comp,ULSnSLS_Group_Design_ten};
%sorted_ULSnSLS_Group_Design = {sorted_ULSnSLS_Group_Design_comp,sorted_ULSnSLS_Group_Design_ten};
number_of_pile_types_ULSnSLS_Group_Design = {pile_number_ULSnSLS_Design_comp, pile_number_ULSnSLS_Design_ten};
ULSnSLS_Group_Cost = {per_pile_cost_ULSnSLS_Design_comp, group_pile_cost_ULSnSLS_Design_comp, per_pile_cost_ULSnSLS_Design_ten, group_pile_cost_ULSnSLS_Design_ten};

disp('ULSnSLS done')
end
