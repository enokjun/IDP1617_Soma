% Copyright (C) 2017 Enok Cheon

% Author: Enok Cheon <Enok Cheon@LAPTOP-I5VE6B9I>
% Created: 2017-03-03

function [dx, xn, yn, Ixx, Iyy, bb, ll, Ab, perimeter, As, N_pile, pile_coordinates, pile_types, which_cell_no] = Y7_ULSnSLS_core_design_N(d,L,N)

%% the purpose of this to retrive the group pile information without over crowding the mother code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input from the ULS analysis - N based on the pile diameter due to space constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assumption: centre-to-centre distance between piles (dx) is 3 x diameter
dx = 2.5*d;       % as 3xd was giving such a big pile cap, try 2.5

% extract coordinates of each piles and its type for each layout
switch N
    case 4 % 2x2
        xn = dx*[-0.5; 0.5; -0.5; 0.5];
        yn = dx*[0.5; 0.5; -0.5; -0.5];
        pile_types = [1; 1; 1; 1];  
        which_cell_no = 1;
    case 6 % 2x3
        xn = dx*[-1; 0; 1; -1; 0; 1];
        yn = dx*[0.5; 0.5; 0.5; -0.5; -0.5; -0.5];
        pile_types = [1; 2; 1; 1; 2; 1];     % 1 = corner, 2 = edge
        which_cell_no = 2;
    case 8 % 2x4
        xn = dx*[-1.5; -0.5; 0.5; 1.5;-1.5; -0.5; 0.5; 1.5];
        yn = dx*[0.5; 0.5; 0.5; 0.5; -0.5; -0.5; -0.5; -0.5];
        pile_types = [1; 2; 2; 1; 1; 2; 2; 1];     % 1 = corner, 2 = edge
        which_cell_no = 3;
    case 9  % 3x3 layout
        xn = dx*[-1; 0; 1; -1; 0; 1; -1; 0; 1];
        yn = dx*[1; 1; 1; 0; 0; 0; -1; -1; -1];
        pile_types = [1; 2; 1; 2; 3; 2; 1; 2; 1];
        which_cell_no = 4;    
    case 10  % 2x5 layout
        xn = dx*[-2; -1; 0; 1; 2; -2; -1; 0; 1; 2];
        yn = dx*[0.5; 0.5; 0.5; 0.5; 0.5; -0.5; -0.5; -0.5; -0.5; -0.5];
        pile_types = [1; 2; 1; 2; 3; 2; 1; 2; 1];   
        which_cell_no = 5;
    case 12 % 3x4
        Layout_3x4 = dlmread('pile_3x4.csv',',',0,0);
        xn = Layout_3x4(:,1)*dx;
        yn = Layout_3x4(:,2)*dx;
        pile_types = Layout_3x4(:,3);
        which_cell_no = 6;
    case 15  % 3x5 layout
        Layout_3x5 = dlmread('pile_3x5.csv',',',0,0);
        xn = Layout_3x5(:,1)*dx;
        yn = Layout_3x5(:,2)*dx;
        pile_types = Layout_3x5(:,3);
        which_cell_no = 7;
    case 16  % 4x4 layout
        Layout_4x4 = dlmread('pile_4x4.csv',',',0,0);
        xn = Layout_4x4(:,1)*dx;
        yn = Layout_4x4(:,2)*dx;
        pile_types = Layout_4x4(:,3);
        which_cell_no = 8;
    case 18  % 3x6 layout
        Layout_3x6 = dlmread('pile_3x6.csv',',',0,0);
        xn = Layout_3x6(:,1)*dx;
        yn = Layout_3x6(:,2)*dx;
        pile_types = Layout_3x6(:,3);
        which_cell_no = 9;
    case 20 % 4x5
        Layout_4x5 = dlmread('pile_4x5.csv',',',0,0);
        xn = Layout_4x5(:,1)*dx;
        yn = Layout_4x5(:,2)*dx;
        pile_types = Layout_4x5(:,3); 
        which_cell_no = 10;
    case 21 % 3x7
        Layout_3x7 = dlmread('pile_3x7.csv',',',0,0);
        xn = Layout_3x7(:,1)*dx;
        yn = Layout_3x7(:,2)*dx;
        pile_types = Layout_3x7(:,3);  
        which_cell_no = 11;
    case 24 % 4x6
        Layout_4x6 = dlmread('pile_4x6.csv',',',0,0);
        xn = Layout_4x6(:,1)*dx;
        yn = Layout_4x6(:,2)*dx;
        pile_types = Layout_4x6(:,3);  
        which_cell_no = 12;
    case 25 % 5x5
        Layout_5x5 = dlmread('pile_5x5.csv',',',0,0);
        xn = Layout_5x5(:,1)*dx;
        yn = Layout_5x5(:,2)*dx;
        pile_types = Layout_5x5(:,3);  
        which_cell_no = 13;
    case 28 % 4x7
        Layout_4x7 = dlmread('pile_4x7.csv',',',0,0);
        xn = Layout_4x7(:,1)*dx;
        yn = Layout_4x7(:,2)*dx;
        pile_types = Layout_4x7(:,3); 
        which_cell_no = 14;          
    case 30 % 5x6
        Layout_5x6 = dlmread('pile_5x6.csv',',',0,0);
        xn = Layout_5x6(:,1)*dx;
        yn = Layout_5x6(:,2)*dx;
        pile_types = Layout_5x6(:,3);   
        which_cell_no = 15;  
    case 32 % 4x8
        Layout_4x8 = dlmread('pile_4x8.csv',',',0,0);
        xn = Layout_4x8(:,1)*dx;
        yn = Layout_4x8(:,2)*dx;
        pile_types = Layout_4x8(:,3);  
        which_cell_no = 16;         
    case 35 % 5x7
        Layout_5x7 = dlmread('pile_5x7.csv',',',0,0);
        xn = Layout_5x7(:,1)*dx;
        yn = Layout_5x7(:,2)*dx;
        pile_types = Layout_5x7(:,3); 
        which_cell_no = 17;    
    case 361 % 4x9
        Layout_4x9 = dlmread('pile_4x9.csv',',',0,0);
        xn = Layout_4x9(:,1)*dx;
        yn = Layout_4x9(:,2)*dx;
        pile_types = Layout_4x9(:,3);  
        which_cell_no = 18;  
    case 362 % 6x6
        Layout_6x6 = dlmread('pile_6x6.csv',',',0,0);
        xn = Layout_6x6(:,1)*dx;
        yn = Layout_6x6(:,2)*dx;
        pile_types = Layout_6x6(:,3);    
        which_cell_no = 19;
    case 40 % 5x8
        Layout_5x8 = dlmread('pile_5x8.csv',',',0,0);
        xn = Layout_5x8(:,1)*dx;
        yn = Layout_5x8(:,2)*dx;
        pile_types = Layout_5x8(:,3);   
        which_cell_no = 20;  
    case 42 % 6x7
        Layout_6x7 = dlmread('pile_6x7.csv',',',0,0);
        xn = Layout_6x7(:,1)*dx;
        yn = Layout_6x7(:,2)*dx;
        pile_types = Layout_6x7(:,3);   
        which_cell_no = 21; 
    case 44 % 4x11
        Layout_4x11 = dlmread('pile_4x11.csv',',',0,0);
        xn = Layout_4x11(:,1)*dx;
        yn = Layout_4x11(:,2)*dx;
        pile_types = Layout_4x11(:,3); 
        which_cell_no = 22;
    case 45 % 5x9
        Layout_5x9 = dlmread('pile_5x9.csv',',',0,0);
        xn = Layout_5x9(:,1)*dx;
        yn = Layout_5x9(:,2)*dx;
        pile_types = Layout_5x9(:,3);   
        which_cell_no = 23;  
    case 48 % 6x8
        Layout_6x8 = dlmread('pile_6x8.csv',',',0,0);
        xn = Layout_6x8(:,1)*dx;
        yn = Layout_6x8(:,2)*dx;
        pile_types = Layout_6x8(:,3);   
        which_cell_no = 24;
    case 49 % 7x7
        Layout_7x7 = dlmread('pile_7x7.csv',',',0,0);
        xn = Layout_7x7(:,1)*dx;
        yn = Layout_7x7(:,2)*dx;
        pile_types = Layout_7x7(:,3);
        which_cell_no = 25; 
    case 50 % 5x10
        Layout_5x10 = dlmread('pile_5x10.csv',',',0,0);
        xn = Layout_5x10(:,1)*dx;
        yn = Layout_5x10(:,2)*dx;
        pile_types = Layout_5x10(:,3);   
        which_cell_no = 26;
    case 54 % 6x9
        Layout_6x9 = dlmread('pile_6x9.csv',',',0,0);
        xn = Layout_6x9(:,1)*dx;
        yn = Layout_6x9(:,2)*dx;
        pile_types = Layout_6x9(:,3); 
        which_cell_no = 27;
    case 55 % 5x11
        Layout_5x11 = dlmread('pile_5x11.csv',',',0,0);
        xn = Layout_5x11(:,1)*dx;
        yn = Layout_5x11(:,2)*dx;
        pile_types = Layout_5x11(:,3); 
        which_cell_no = 28;
    case 56 % 7x8
        Layout_7x8 = dlmread('pile_7x8.csv',',',0,0);
        xn = Layout_7x8(:,1)*dx;
        yn = Layout_7x8(:,2)*dx;
        pile_types = Layout_7x8(:,3); 
        which_cell_no = 29;
    case 60 % 6x10
        Layout_6x10 = dlmread('pile_6x10.csv',',',0,0);
        xn = Layout_6x10(:,1)*dx;
        yn = Layout_6x10(:,2)*dx;
        pile_types = Layout_6x10(:,3); 
        which_cell_no = 30;
    case 63 % 7x9
        Layout_7x9 = dlmread('pile_7x9.csv',',',0,0);
        xn = Layout_7x9(:,1)*dx;
        yn = Layout_7x9(:,2)*dx;
        pile_types = Layout_7x9(:,3); 
        which_cell_no = 31;
    case 64 % 8x8
        Layout_8x8 = dlmread('pile_8x8.csv',',',0,0);
        xn = Layout_8x8(:,1)*dx;
        yn = Layout_8x8(:,2)*dx;
        pile_types = Layout_8x8(:,3); 
        which_cell_no = 32;
    case 66 % 6x11
        Layout_6x11 = dlmread('pile_6x11.csv',',',0,0);
        xn = Layout_6x11(:,1)*dx;
        yn = Layout_6x11(:,2)*dx;
        pile_types = Layout_6x11(:,3); 
        which_cell_no = 33;
    case 70 % 7x10
        Layout_7x10 = dlmread('pile_7x10.csv',',',0,0);
        xn = Layout_7x10(:,1)*dx;
        yn = Layout_7x10(:,2)*dx;
        pile_types = Layout_7x10(:,3); 
        which_cell_no = 34;
    case 72 % 8x9
        Layout_8x9 = dlmread('pile_8x9.csv',',',0,0);
        xn = Layout_8x9(:,1)*dx;
        yn = Layout_8x9(:,2)*dx;
        pile_types = Layout_8x9(:,3);
        which_cell_no = 35;
    case 77 % 7x11
        Layout_7x11 = dlmread('pile_7x11.csv',',',0,0);
        xn = Layout_7x11(:,1)*dx;
        yn = Layout_7x11(:,2)*dx;
        pile_types = Layout_7x11(:,3); 
        which_cell_no = 36;
    case 80 % 8x10
        Layout_8x10 = dlmread('pile_8x10.csv',',',0,0);
        xn = Layout_8x10(:,1)*dx;
        yn = Layout_8x10(:,2)*dx;
        pile_types = Layout_8x10(:,3); 
        which_cell_no = 37;
    case 88 % 8x11
        Layout_8x11 = dlmread('pile_8x11.csv',',',0,0);
        xn = Layout_8x11(:,1)*dx;
        yn = Layout_8x11(:,2)*dx;
        pile_types = Layout_8x11(:,3);
        which_cell_no = 38; 
    otherwise
        disp(N)
        error('check your codes or input at line pile group selection part')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate design parameters from coordinates and other values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for_bb = min(numberofpilesinx, numberofpilesiny);
for_ll = max(numberofpilesinx, numberofpilesiny);
bb = d+(for_bb-1)*dx;
ll = d+(for_ll-1)*dx;
Ab = bb*ll;      
perimeter = 2*(bb+ll);
As = L*perimeter;

%% distribution of axial loads
pile_coordinates = [xn, yn];       % coordinate showing pile locations
N_pile = [];
unique_piles = unique(pile_types);
for unique_piles_analyse = 1:length(unique_piles)
    N_pile(unique_piles_analyse) = [sum(pile_types == unique_piles(unique_piles_analyse))];
end

end
