% Copyright (C) 2016 Enok Cheon
% Author: Enok Cheon <Enok Cheon@LAPTOP-I5VE6B9I>
% Created: 2016-12-31

function [sorted_ULSnSLS_Design, ULSnSLS_Design] = Z5_ULSnSLS_CFA_wall(P_list, ciria_data, STR_GEO_AMR, Beta_data, conc_data, soil, basement_data, pile_dim, sett_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads
Dead_P = P_list(1);
Live_P = P_list(2);
Wind_P = P_list(3);

% pile dimension limits
D = pile_dim{1};
Lmin = pile_dim{2};
Lmax = pile_dim{3};
Lint = pile_dim{4};

% basement data
basement_depth = basement_data(2);
soil_otherside_depth = basement_data(4);

% beta method
Ks = Beta_data(1);
del_fill = Beta_data(2);
del_tergravel = Beta_data(3);

% concrete data
fck_cube = conc_data(1);
Ep = conc_data(2);

% CIRIA data
Mmin = ciria_data{1};
Mmax = ciria_data{2};
Space = ciria_data{3};

% settlement data
Ms = sett_data(1);
Ke = sett_data(2);
Lo = sett_data(3);
slim_single = sett_data(4);
slim_group = sett_data(5);
slim_percent_single = sett_data(6);
slim_percent_group = sett_data(7); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLS calculation iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
comb = {[1, 3, 5],[2, 4, 6]};       % factors to extract from STR_GEO_AMR. (1) is Comb 1 and (2) is Comb 2

% matrix to store analysed values
ULSnSLS_design_C1_comp = [];
ULSnSLS_design_C2_comp = [];
ULSnSLS_design_C1_ten = [];
ULSnSLS_design_C2_ten = [];
ULSnSLS_Design_comp = [];
ULSnSLS_Design_ten = [];
sorted_ULSnSLS_Design_comp = [];
sorted_ULSnSLS_Design_ten = [];

results = [];
tab = 1;


for i = 1:2
    row = 1;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load combination analysis - total axial loads
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% factored loads
    ULS_Dead_Unf = STR_GEO_AMR{comb{i}(1)}(1)*Dead_P;
    ULS_Dead_F = STR_GEO_AMR{comb{i}(1)}(2)*Dead_P;
    ULS_Live_Unf = STR_GEO_AMR{comb{i}(1)}(3)*Live_P;
    ULS_Wind_Unf = STR_GEO_AMR{comb{i}(1)}(4)*Wind_P;
    SLS_Dead_Unf = Dead_P;
    SLS_Dead_F = Dead_P;
    SLS_Live_Unf = Live_P;
    SLS_Wind_Unf = Wind_P;

    ULS_Load_comp_total = 0;
    ULS_Load_ten_total = 0;
    SLS_Load_comp_total = 0;
    SLS_Load_ten_total = 0;

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
    end
    %% dead and live load (no wind)
    if (Live_P == 0 && Dead_P > 0)                                       % no live
        ULS_load_2_comp_total = 0;
        ULS_load_2_ten_total = 0;
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
    % combination case 1, 3
    elseif (Live_P  > 0 && Dead_P < 0) || (Live_P < 0 && Dead_P > 0)         % gravity only - different load type for each loads
        % axial
        ULS_load_2a = ULS_Dead_Unf(1);                          % live is favourable; therefore, should not be considered
        ULS_load_2b = ULS_Dead_F(1) + ULS_Live_Unf(1);          % dead is favourable
        ULS_load_2_comp_total = max(ULS_load_2a,ULS_load_2b);
        ULS_load_2_ten_total = min(ULS_load_2a,ULS_load_2b);
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
    % combination case 5,1
    elseif (Live_P> 0 && Dead_P < 0 && Wind_P > 0)||(Live_P< 0 && Dead_P > 0 && Wind_P < 0)   % live and wind are same; dead is different 
        ULS_load_3a = ULS_Dead_F(1) + ULS_Live_Unf(1) + ULS_Wind_Unf(1);
        ULS_load_3b = ULS_Dead_Unf(1);
        ULS_load_3_comp_total = max(ULS_load_3a,ULS_load_3b);
        ULS_load_3_ten_total = min(ULS_load_3a,ULS_load_3b);

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
       
    % combination case 6,3
    elseif (Live_P < 0 && Dead_P > 0 && Wind_P > 0)||(Live_P> 0 && Dead_P < 0 && Wind_P < 0)   % dead and wind are same; live is different
        ULS_load_3a = ULS_Dead_Unf(1) + ULS_Wind_Unf(1);
        ULS_load_3b = ULS_Dead_F(1) + ULS_Live_Unf(1);  
        ULS_load_3_comp_total = max(ULS_load_3a,ULS_load_3b);
        ULS_load_3_ten_total = min(ULS_load_3a,ULS_load_3b);

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
        case 91
            SLS_Load_comp_total = 0;
            SLS_Load_ten_total = SLS_Dead_Unf(1); 
        case 13
            SLS_Load_comp_total = SLS_Dead_Unf(1); 
            SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
        case 31
            SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
            SLS_Load_ten_total = SLS_Dead_Unf(1); 
        case 15
            SLS_Load_comp_total = SLS_Dead_Unf(1); 
            SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
        case 51
            SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1); 
            SLS_Load_ten_total = SLS_Dead_Unf(1); 
        case 29
            SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
            SLS_Load_ten_total = 0;
        case 92
            SLS_Load_comp_total = 0;
            SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1);
        case 49
            SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
            SLS_Load_ten_total = 0;
        case 94
            SLS_Load_comp_total = 0;
            SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Live_Unf(1)+SLS_Wind_Unf(1);
        case 36
            SLS_Load_comp_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
            SLS_Load_ten_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
        case 63
            SLS_Load_comp_total = SLS_Dead_Unf(1)+SLS_Wind_Unf(1);
            SLS_Load_ten_total = SLS_Dead_F(1)+SLS_Live_Unf(1); 
    end

    % list of total axial loads to design
    if isempty(SLS_Load_comp_total) == 1 || isempty(SLS_Load_ten_total) == 1 % no loads
        error('something wrong with the loads')
    elseif (isempty(SLS_Load_comp_total) == 0 && isempty(SLS_Load_ten_total) == 0)  % both tension and compression occurs
        SLS_Load_total_list = [SLS_Load_comp_total, SLS_Load_ten_total];
    end

    for d = D
        % determine space in between
        if d == 0.6
            M_space = Space(1);
        elseif d == 0.75
            M_space = Space(2);
        elseif d == 0.9
            M_space = Space(3);
        elseif d == 1.05
            M_space = Space(4);
        elseif d == 1.2
            M_space = Space(5);
        else
            error('Something wrong with diameter');
        end

        for Z = Lmin: Lint: Lmax
            for M = Mmin: Mmax
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% iterate SLS_Group_Design
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for loads_cases = 1:2
                    ULS_load_total = ULS_Load_total_list(loads_cases);
                    SLS_load_total = SLS_Load_total_list(loads_cases);
                    if SLS_load_total == 0  || ULS_load_total == 0 
                        if loads_cases == 1
                            continue
                        elseif loads_cases == 2
                            break
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    %% base resistance stress
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % base on terrace gravel
                    qb = 0;
                    if ULSnSLS_soil_Lonclay1(1) > Z && Z >= ULSnSLS_soil_tergravel(1)
                        Nq = exp(pi*tand(ULSnSLS_soil_tergravel(5)))*(tand(45+0.5*ULSnSLS_soil_tergravel(5))).^2;
                        Nq_SLS = exp(pi*tand(soil(3,5)))*(tand(45+0.5*soil(3,5))).^2;
                        % overburden stress
                        overburden_stress = (d-ULSnSLS_soil_tergravel(1))*ULSnSLS_soil_tergravel(3)+ ...
                            ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+ (ULSnSLS_soil_fill(1)-basement_depth)*ULSnSLS_soil_fill(3);
                        overburden_stress_SLS = (d-soil(3,1))*soil(3,3)+soil(2,2)*soil(2,3)+(soil(1,1)-basement_depth)*soil(1,3);
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
                        qb_SLS = 9*soil(6,4);
                    end
                    % base
                    wall_length = d+(M-1)*M_space;
                    Ab = d*wall_length;         % shaft cross-sectional area
                    ULS_Rb = Ab*qb;
                    SLS_Rb = Ab*qb_SLS;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% shaft resistance stress
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % shaft resistance depth and area
                    t_shaft_dep = [basement_depth, Z; soil_otherside_depth, Z];
                    len_fill = [];
                    len_orgclay = [];
                    len_tergravel = [];
                    len_Lonclay1 = [];
                    len_Lonclay2 = [];
                    len_Lonclay3 = [];
                    for length_of_piles = 1:2
                        %% shaft on fill
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_fill(1)
                            len_fill(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_fill(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_fill(1) && t_shaft_dep(length_of_piles,2) < ULSnSLS_soil_orgclay(1)
                                len_fill(length_of_piles) = t_shaft_dep(length_of_piles,2);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_fill(1) && t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_orgclay(1)
                                len_fill(length_of_piles) = ULSnSLS_soil_fill(1) + ULSnSLS_soil_fill(2);
                            else
                                len_fill(length_of_piles) = 0;
                            end
                        end
                        % beta
                        beta_fill = Ks*tand(del_fill);
                        % pressure
                        qs_fill = beta_fill*(len_fill(length_of_piles)*0.5)*ULSnSLS_soil_fill(3);
                        qs_fill_SLS = beta_fill*(len_fill(length_of_piles)*0.5)*soil(1,3);
                        
                        %% shaft on organic clay
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_orgclay(1)
                            len_orgclay(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_orgclay(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_orgclay(1) && t_shaft_dep(length_of_piles,2) < ULSnSLS_soil_tergravel(1)
                                len_orgclay(length_of_piles) = t_shaft_dep(length_of_piles,2)-ULSnSLS_soil_orgclay(1);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_orgclay(1) && t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_tergravel(1)
                                len_orgclay(length_of_piles) = ULSnSLS_soil_orgclay(2);
                            else
                                len_orgclay(length_of_piles) = 0;
                            end  
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
                        
                        %% shaft on terrace gravel
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_tergravel(1)
                            len_tergravel(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_tergravel(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_tergravel(1) && t_shaft_dep(length_of_piles,2) < ULSnSLS_soil_Lonclay1(1)
                                len_tergravel(length_of_piles) = t_shaft_dep(length_of_piles,2)-ULSnSLS_soil_tergravel(1);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_tergravel(1) && t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay1(1)
                                len_tergravel(length_of_piles) = ULSnSLS_soil_tergravel(2);
                            else
                                len_tergravel(length_of_piles) = 0;
                            end    
                        end
                        % beta
                        beta_tergravel = Ks*tand(del_tergravel);
                        % pressure
                        overburden_pres = (len_tergravel(length_of_piles)*0.5)*ULSnSLS_soil_tergravel(3)+...
                                           ULSnSLS_soil_orgclay(2)*ULSnSLS_soil_orgclay(3)+...
                                          (ULSnSLS_soil_fill(1)+ULSnSLS_soil_fill(2)-t_shaft_dep(length_of_piles,1))*ULSnSLS_soil_fill(3);
                        qs_tergravel = beta_tergravel*overburden_pres;
                        overburden_pres_SLS = (len_tergravel(length_of_piles)*0.5)*soil(3,3)+soil(2,2)*soil(2,3)+(soil(1,1)+soil(1,2)-t_shaft_dep(length_of_piles,1))*soil(1,3);
                        qs_tergravel_SLS = beta_tergravel*overburden_pres_SLS;
                        
                        %% shaft on London clay
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_Lonclay1(1)
                            len_Lonclay1(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_Lonclay1(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep(length_of_piles,2) < ULSnSLS_soil_Lonclay2(1)
                                len_Lonclay1(length_of_piles) = t_shaft_dep(length_of_piles,2)-ULSnSLS_soil_Lonclay1(1);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay1(1) && t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay2(1)
                                len_Lonclay1(length_of_piles) = ULSnSLS_soil_Lonclay1(2);
                            else
                                len_Lonclay1(length_of_piles) = 0;
                            end   
                        end 
                        % pressure
                        if ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4) > 110
                            qs_Lonclay1 = 110;
                        else
                            qs_Lonclay1 = ULSnSLS_soil_Lonclay1(7)*ULSnSLS_soil_Lonclay1(4);
                        end
                        if soil(4,7)*soil(4,4) > 110
                            qs_Lonclay1_SLS = 110;
                        else
                            qs_Lonclay1_SLS = soil(4,7)*soil(4,4) ;
                        end
                        
                        %% shaft on Stiff London clay
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_Lonclay2(1)
                            len_Lonclay2(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_Lonclay2(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep(length_of_piles,2) < ULSnSLS_soil_Lonclay3(1)
                                len_Lonclay2(length_of_piles) = t_shaft_dep(length_of_piles,2)-ULSnSLS_soil_Lonclay2(1);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay2(1) && t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay3(1)
                                len_Lonclay2(length_of_piles) = ULSnSLS_soil_Lonclay2(2);
                            else
                                len_Lonclay2(length_of_piles) = 0;
                            end   
                        end 
                        % pressure
                        if ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4) > 110
                            qs_Lonclay2 = 110;
                        else
                            qs_Lonclay2 = ULSnSLS_soil_Lonclay2(7)*ULSnSLS_soil_Lonclay2(4);
                        end
                        if soil(5,7)*soil(5,4) > 110
                            qs_Lonclay2_SLS = 110;
                        else
                            qs_Lonclay2_SLS = soil(5,7)*soil(5,4);
                        end
                        
                        %% shaft on Claystone
                        % length
                        if t_shaft_dep(length_of_piles,1) >= ULSnSLS_soil_Lonclay3(1)
                            len_Lonclay3(length_of_piles) = 0;
                        elseif t_shaft_dep(length_of_piles,1) < ULSnSLS_soil_Lonclay3(1)
                            if t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep(length_of_piles,2) < (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
                                len_Lonclay3(length_of_piles) = t_shaft_dep(length_of_piles,2)-ULSnSLS_soil_Lonclay3(1);
                            elseif t_shaft_dep(length_of_piles,2) >= ULSnSLS_soil_Lonclay3(1) && t_shaft_dep(length_of_piles,2) >= (ULSnSLS_soil_Lonclay3(1)+ULSnSLS_soil_Lonclay3(2))
                                len_Lonclay3(length_of_piles) = ULSnSLS_soil_Lonclay3(2);
                            else
                                len_Lonclay3(length_of_piles) = 0;
                            end 
                        end
                        % pressure
                        if ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4) > 110
                            qs_Lonclay3 = 110;
                        else
                            qs_Lonclay3 = ULSnSLS_soil_Lonclay3(7)*ULSnSLS_soil_Lonclay3(4);
                        end
                        if soil(6,7)*soil(6,4) > 110
                            qs_Lonclay3_SLS = 110;
                        else
                            qs_Lonclay3_SLS = soil(6,7)*soil(6,4);
                        end
                    end

                    if ULS_load_total > 0 
                        % shaft
                        Qs_fill = qs_fill*wall_length*sum(len_fill);
                        Qs_org_clay = qs_org_clay*wall_length*sum(len_orgclay);
                        Qs_tergravel = qs_tergravel*wall_length*sum(len_tergravel);
                        Qs_Lonclay1 = qs_Lonclay1*wall_length*sum(len_Lonclay1);
                        Qs_Lonclay2 = qs_Lonclay2*wall_length*sum(len_Lonclay2);
                        Qs_Lonclay3 = qs_Lonclay3*wall_length*sum(len_Lonclay3);
                        ULS_Rs = -Qs_fill-Qs_org_clay+Qs_tergravel+Qs_Lonclay1+Qs_Lonclay2+Qs_Lonclay3;
                        % negative skin friction at fill and organic clay layers

                        Qs_fill_SLS = qs_fill_SLS*wall_length*sum(len_fill);
                        Qs_org_clay_SLS = qs_org_clay_SLS*wall_length*sum(len_orgclay);
                        Qs_tergravel_SLS = qs_tergravel_SLS*wall_length*sum(len_tergravel);
                        Qs_Lonclay1_SLS = qs_Lonclay1_SLS*wall_length*sum(len_Lonclay1);
                        Qs_Lonclay2_SLS = qs_Lonclay2_SLS*wall_length*sum(len_Lonclay2);
                        Qs_Lonclay3_SLS = qs_Lonclay3_SLS*wall_length*sum(len_Lonclay3);
                        SLS_Rs = -Qs_fill_SLS-Qs_org_clay_SLS+Qs_tergravel_SLS+Qs_Lonclay1_SLS+Qs_Lonclay2_SLS+Qs_Lonclay3_SLS;
                        % negative skin friction at fill and organic clay layers

                        % Ultimate resistance and check (LSDA Foundation No 1, 3.2) - compression
                        ULS_Ru_from_Rs_Rb = ULS_Rb/STR_GEO_AMR{comb{i}(3)}(1) + ULS_Rs/STR_GEO_AMR{comb{i}(3)}(2);    % unit: kN
                        ULS_Ru_from_Rs = ULS_Rs/1.2;                               % unit: kN
                        ULS_Ru_from_fck_cube = Ab*fck_cube*0.25;                    % unit: kN
                        ULS_Ru_comp = min([ULS_Ru_from_Rs_Rb;ULS_Ru_from_Rs;ULS_Ru_from_fck_cube]);
                        SLS_Ru_comp = SLS_Rb+SLS_Rs;

                    elseif ULS_load_total < 0                 
                        % shaft
                        % assume arc length of 55% of the circumference
                        Qs_fill = qs_fill*(0.55*pi*d)*sum(len_fill);
                        Qs_org_clay = qs_org_clay*(0.55*pi*d)*sum(len_orgclay);
                        Qs_tergravel = qs_tergravel*(0.55*pi*d)*sum(len_tergravel);
                        Qs_Lonclay1 = qs_Lonclay1*(0.55*pi*d)*sum(len_Lonclay1);
                        Qs_Lonclay2 = qs_Lonclay2*(0.55*pi*d)*sum(len_Lonclay2);
                        Qs_Lonclay3 = qs_Lonclay3*(0.55*pi*d)*sum(len_Lonclay3);
                        ULS_Rs = -Qs_fill-Qs_org_clay+Qs_tergravel+Qs_Lonclay1+Qs_Lonclay2+Qs_Lonclay3;
                        % negative skin friction at fill and organic clay layers

                        Qs_fill_SLS = qs_fill_SLS*(0.55*pi*d)*sum(len_fill);
                        Qs_org_clay_SLS = qs_org_clay_SLS*(0.55*pi*d)*sum(len_orgclay);
                        Qs_tergravel_SLS = qs_tergravel_SLS*(0.55*pi*d)*sum(len_tergravel);
                        Qs_Lonclay1_SLS = qs_Lonclay1_SLS*(0.55*pi*d)*sum(len_Lonclay1);
                        Qs_Lonclay2_SLS = qs_Lonclay2_SLS*(0.55*pi*d)*sum(len_Lonclay2);
                        Qs_Lonclay3_SLS = qs_Lonclay3_SLS*(0.55*pi*d)*sum(len_Lonclay3);
                        SLS_Rs = -Qs_fill_SLS-Qs_org_clay_SLS+Qs_tergravel_SLS+Qs_Lonclay1_SLS+Qs_Lonclay2_SLS+Qs_Lonclay3_SLS;
                        % negative skin friction at fill and organic clay layers
                    
                        % Ultimate resistance -  tension
                        ULS_Ru_ten = M*ULS_Rs/STR_GEO_AMR{comb{i}(3)}(3);
                        SLS_Ru_ten = M*SLS_Rs;
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Settlement   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                    Q = SLS_load_total;
                    if SLS_load_total > 0
                        SLS_Ru = SLS_Ru_comp;
                    elseif SLS_load_total < 0  
                        SLS_Ru = SLS_Ru_ten;
                    end

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
                    alpha = SLS_Rs;
                    if Q > 0
                        beta = SLS_Rb*d*Eb;
                        delta = SLS_Rb*0.6;
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
                    if abs(Q) <= abs(SLS_Rs)
                        se = 1000*(4*Q*(Lo + Ke*(Z-Lo)))/(pi*Ep*(d.^2));
                    elseif abs(Q) > abs(SLS_Rs) && abs(Q) <= abs(SLS_Ru)
                        se = 1000*(4*((Q*Z)-(SLS_Rs*(1-Ke)*(Z-Lo))))/(pi*Ep*(d.^2));
                    elseif abs(Q) > abs(SLS_Ru)
                        se = 9999;
                    end
                    
                    % modify datalog
                    if se==9999 || sr==9999
                        s = 9999;
                    else
                        s = sr + se;
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %% Record and tabulate
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    % compression
                    if (ULS_load_total > 0 && ULS_load_total <= ULS_Ru_comp) && (abs(s) <= (slim_percent_single/100)*d*1000)
                        % Material Volume
                        V = d*wall_length*Z;
                        N = 2*M-1;

                        % record
                        dimensions = [Z, d, wall_length, M, N]; 
                        loads = [critical_load_combi, ULS_load_total, SLS_load_total]; 
                        resistance = [ULS_Rb, ULS_Rs, ULS_Ru_from_Rs, ULS_Ru_from_fck_cube, ULS_Ru_comp, SLS_Rb, SLS_Rs, SLS_Ru_comp]; 
                        settlements = [sr, se, s]; 
                        stiffness = [SLS_load_total/(0.001*s)];
                        costs = [V];
                        
                        if i == 1
                            ULSnSLS_design_C1_comp(row,:) = [dimensions, loads, resistance, settlements, stiffness, costs];
                            row = row+1;
                        elseif i == 2
                            ULSnSLS_design_C2_comp(row,:) = [dimensions, loads, resistance, settlements, stiffness, costs];
                            row = row+1;
                        end

                        results(tab,:) = [dimensions, loads, resistance, settlements, costs];
                        tab = tab+1;

                    % tension
                    elseif (ULS_load_total < 0 && abs(ULS_load_total) <= abs(ULS_Ru_ten)) && (abs(s) <= (slim_percent_single/100)*d*1000)
                        % Material Volume
                        V = d*wall_length*Z;
                        N = 2*M-1;

                        % record
                        dimensions = [Z, d, wall_length, M, N]; 
                        loads = [critical_load_combi, ULS_load_total, SLS_load_total]; 
                        resistance = [ULS_Rs, -ULS_Ru_ten, SLS_Rs, -SLS_Ru_ten]; 
                        settlements = [sr, se, s]; 
                        stiffness = [SLS_load_total*1000/s];  % unit: kN/m
                        costs = [V];
                        
                        if i == 1
                            ULSnSLS_design_C1_ten(row,:) = [dimensions, loads, resistance, settlements, stiffness, costs];
                            row = row+1;
                        elseif i == 2
                            ULSnSLS_design_C2_ten(row,:) = [dimensions, loads, resistance, settlements, stiffness, costs];
                            row = row+1;
                        end

                        results(tab,:) = [dimensions, loads, resistance, settlements, stiffness, costs];
                        tab = tab+1;
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
    %% compare load results
    if ULS_Load_total_list(loads_cases) > 0 && (isempty(ULSnSLS_design_C1_comp) == 0 && isempty(ULSnSLS_design_C2_comp) == 0)
        [w1,~]=size(ULSnSLS_design_C1_comp);
        [w2,~]=size(ULSnSLS_design_C2_comp);
        % these variables stores the number of times a combination gives a worse scenario
        for l = 1:w1
            for m = n:w2
                if ULSnSLS_design_C1_comp(l,1) == ULSnSLS_design_C2_comp(m,1) && ULSnSLS_design_C1_comp(l,2) == ULSnSLS_design_C2_comp(m,2) && ULSnSLS_design_C1_comp(l,4) == ULSnSLS_design_C2_comp(m,4)
                    if ULSnSLS_design_C1_comp(l,13) >= ULSnSLS_design_C2_comp(m,13) 
                        ULSnSLS_Design_comp(kk,:) = [ULSnSLS_design_C2_comp(m,:), 2]; 
                        kk = kk + 1;
                        n = m;
                    elseif ULSnSLS_design_C1_comp(l,13) < ULSnSLS_design_C2_comp(m,13) 
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
                if ULSnSLS_design_C1_ten(l,1) == ULSnSLS_design_C2_ten(m,1) && ULSnSLS_design_C1_ten(l,2) == ULSnSLS_design_C2_ten(m,2) && ULSnSLS_design_C1_ten(l,4) == ULSnSLS_design_C2_ten(m,4)
                    if abs(ULSnSLS_design_C1_ten(l,10)) >= abs(ULSnSLS_design_C2_ten(m,10))
                        ULSnSLS_Design_ten(kk,:) = [ULSnSLS_design_C2_ten(m,:), 2]; 
                        kk = kk + 1;
                        n = m;
                    elseif abs(ULSnSLS_design_C1_ten(l,13)) < abs(ULSnSLS_design_C2_ten(m,13)) 
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

ULSnSLS_Design = {ULSnSLS_design_C1_comp, ULSnSLS_design_C2_comp, ULSnSLS_design_C1_ten, ULSnSLS_design_C2_ten, ULSnSLS_Design_comp, ULSnSLS_Design_ten};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort accepted design according to Volume of Concrete used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ULSnSLS_Design{5}(1,1) == 9999
    sorted_ULSnSLS_Design_comp = [9999];
else 
    sorted_ULSnSLS_Design_comp = ULSnSLS_Design{5};
    [~,order] = sort(sorted_ULSnSLS_Design_comp(:,end-1));
    sorted_ULSnSLS_Design_comp = sorted_ULSnSLS_Design_comp(order,:);
end

if ULSnSLS_Design{6}(1,1) == 9999
    sorted_ULSnSLS_Design_ten = [9999];
else 
    sorted_ULSnSLS_Design_ten = ULSnSLS_Design{6};
    [~,order] = sort(sorted_ULSnSLS_Design_ten(:,end-1));
    sorted_ULSnSLS_Design_ten = sorted_ULSnSLS_Design_ten(order,:);
end

sorted_ULSnSLS_Design = {sorted_ULSnSLS_Design_comp, sorted_ULSnSLS_Design_ten};

end
