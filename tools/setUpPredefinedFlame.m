function [ p ] = setUpPredefinedFlame( flameName )
%SETUPPREDEFINEDFLAME Sets up a predefined name which is selescted by its
%  name
%
%   Flame available:
%         Cuquel_Cr44_3D
%         Cuquel_Cr44_2D
%         Cuquel_Cr81_3D
%         Cuquel_Cr81_2D       
%         Kornilov_2D
%         Kornilov_3D
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////


switch flameName
  case 'Cuquel_Cr44_3D'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 25e-3;
    s_l_u = 0.27819;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6.7;
    
  case 'Cuquel_Cr44_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 25e-3;
    s_l_u = 0.27819;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6.7;
  
  case 'Cuquel_Cr60_3D'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 18.33e-3;
    s_l_u = 0.27819;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6.7;
    
  case 'Cuquel_Cr60_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 18.33e-3;
    s_l_u = 0.27819;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6.7;
    
  case 'Cuquel_Cr275_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 40e-3;
    s_l_u = 0.37;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 7.155;
    
  case 'Cuquel_Cr81_3D'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 13.55e-3;
    s_l_u = 0.37;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 7.155;
    
  case 'Cuquel_Cr81_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 11e-3;
    R_a = 13.55e-3;
    s_l_u = 0.37;
    u_1_bulkFeed = 0.9;
    u_1_centerIn = u_1_bulkFeed;
    E  = 7.155;
    
  case 'Kornilov_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 1e-3;
    R_a = 2.5e-3;
    s_l_u = 0.37;
    u_1_bulkFeed = 1;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6;
    
  case 'Kornilov_3D'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 1e-3;
    R_a = 2.5e-3;
    s_l_u = 0.37;
    u_1_bulkFeed = 1;
    u_1_centerIn = u_1_bulkFeed;
    E  = 6;
    
  case 'Steinbacher_3D'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 5e-3;
    R_a = 12.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.11;
    E  = 1996/300;
  
  case 'Steinbacher_3D_Cr66'
    Fdim = '3D';
    type = 'inverseV';
    R_i = 5e-3;
    R_a = 7.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.125;
    E  = 1996/300;
    
  case 'Steinbacher_2D'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 5e-3;
    R_a = 12.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.04;
    E  = 1996/300;
  
  case 'Steinbacher_2D_Cr66'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 5e-3;
    R_a = 7.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.07;
    E  = 1996/300;
  
  case 'Steinbacher_2D_Cr17'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 5e-3;
    R_a = 29e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1;
    E  = 1996/300;
    
  case 'Steinbacher_3D_V'
    Fdim = '3D';
    type = 'V';
    R_i = 1.5e-3;
    R_a = 6.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.057;
    E  = 1996/300;
  
  case 'Steinbacher_3D_V2'
    Fdim = '3D';
    type = 'V';
    R_i = 3.5e-3;
    R_a = 8.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = 1.052;
    E  = 1996/300;
  
  case 'Steinbacher_3D_V3'
    Fdim = '3D';
    type = 'V';
    R_i = 5.5e-3;
    R_a = 10.5e-3;
    s_l_u = 0.278;
    u_1_bulkFeed = 1;
    u_1_centerIn = u_1_bulkFeed;
    E  = 1996/300;
    
  case 'SantoshSujith05_2D_u1'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 25e-3;
    C_r = 0.5;
    R_a = R_i/C_r;
    s_l_u = 0.2;
    u_1_bulkFeed = 2.21 * s_l_u;
    u_1_centerIn = u_1_bulkFeed;
    E  = 1000/300;
    
  case 'SantoshSujith05_2D_u2'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 25e-3;
    C_r = 0.5;
    R_a = R_i/C_r;
    s_l_u = 0.2;
    u_1_bulkFeed = 2.98 * s_l_u;
    u_1_centerIn = u_1_bulkFeed;
    E  = 1000/300;
    
  case 'SantoshSujith05_2D_u3'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 25e-3;
    C_r = 0.5;
    R_a = R_i/C_r;
    s_l_u = 0.2;
    u_1_bulkFeed = 4.12 * s_l_u;
    u_1_centerIn = u_1_bulkFeed;
    E  = 1000/300;
    
  case 'Lieuwen99'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 0.1;
    R_a = 0.2;
    s_l_u = 0.278;
    u_1_bulkFeed = 0.55;
    u_1_centerIn = u_1_bulkFeed;
    E  = 2000/300;
    
  case 'Albayrak_V'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 17e-3;
    R_a = 100e-3;
    u_1_bulkFeed = 1;
    s_l_u = sin(pi/4)*u_1_bulkFeed;
    u_1_centerIn = u_1_bulkFeed;
    E  = 1996/300;
  
  case 'Axel_DuctFlame'
    Fdim = '2D';
    type = 'inverseV';
    R_i = 0.1;
    R_a = 0.1;
    s_l_u = 0.278;
    u_1_bulkFeed = 0.55;
    u_1_centerIn = u_1_bulkFeed;
    E  = 2000/300;
    
  otherwise
    error('Unknown flame!')
end

% Set parameters for selected flame
[ p ] = setROMParameters( Fdim , type , R_i , R_a , s_l_u , u_1_bulkFeed , 'Name' , flameName, 'E' , E , ...
  'u1CenterIn' , u_1_centerIn );

end

