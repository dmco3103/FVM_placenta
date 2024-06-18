
%% Finite volume implementation of diffusion in a 3D structured mesh grid %%
% NEEDS GIBBON TOOLBOX FOR PLOTTING, if no Gibbon is there, the plots can be muted

%clc; clear all;

% DOUBLE CHECK UNITS ACROSS SCRIPT!!!!!!!

% PARAMETERS
% distance, surface area, flow rate, velocity: mm, mm^2, mm^3/s, mm/s
% blood viscosity: 3e-3 Pa*s
% permeability: 10^-8 m2 -- 0.01 mm2

% inlet blood pressure: 
inlet_P = 490; % Pin = 2506 Pa (term placenta); = 490 Pa (20weeks)
% outlet blood pressure: 
outlet_P = 0.147;

load 20weeks_cluster.mat

%% DEFINE type of solving we want: stationary or time-dependent
solving = 'stationary';

%% DEFINE type of law to solve
%law = 'diffusion';
law = 'darcy';


%% set up grid

% the mesh grid is defined based on the size of the network being uploaded
% here, a sample terminal fetal tree is used .... TO BE TESTED WITH OTHER
% GEOMETRIES

%load sample_terminal_tree.mat

nb_cell_length=30; %x
nb_cell_height=30; %y
nb_cell_width=30; %z

nb_cell=[nb_cell_length nb_cell_height nb_cell_width]; % number of cells in each direction

  tissue_x=max(tree.X)-min(tree.X)+2*max(tree.D);
  tissue_y=max(tree.Y)-min(tree.Y)+2*max(tree.D);
  tissue_z=max(tree.Z)-min(tree.Z)+2*max(tree.D);

%tissue_x = 40;
%tissue_y = 40;
%tissue_z = 40;

tissue_dimensions=[tissue_x tissue_y tissue_z];

network.center_position=[(min(tree.X)+max(tree.X))/2 (min(tree.Y)+max(tree.Y))/2 (min(tree.Z)+max(tree.Z))/2];

%network.center_position=[0 0 0];

grid=generateGrid2(network.center_position,tissue_dimensions, nb_cell);


%% create mesh with GIBBON function (purely for plotting purposes)
% -- mute this bit if no Gibbon installed

boxDim=[tissue_x tissue_y tissue_z];
%boxDim=[12.7 12.7 24];

boxEl=[nb_cell_length nb_cell_height nb_cell_width];
[meshStruct]=hexMeshBox(boxDim,boxEl,2);
meshStruct.nodes(:,1) = meshStruct.nodes(:,1) + network.center_position(1);
meshStruct.nodes(:,2) = meshStruct.nodes(:,2) + network.center_position(2);
meshStruct.nodes(:,3) = meshStruct.nodes(:,3) + network.center_position(3);


%% Set up variables (add needed information for definition of diffusion matrix later)

% define K permeability tensor (Kxx, Kyy, Kzz)
K = [1e-4 0 0; 0 1e-4 0; 0 0 1e-4]; %same permeability (mm^2) in all directions for now
%K = [0.01 0 0; 0 0.01 0; 0 0 0.01]; %same permeability in all directions for now


switch law

case 'diffusion'
    disp('diffusion')

% add surface areas of each cell face and space steps between centroids

for i=1:length(grid.cell)

%for j=1:length(grid.cell(i).faces)

grid.cell(i).faces.D = {1,1,1,1,1,1};  

grid.cell(i).faces.A = {grid.cell_surfaceXZ,grid.cell_surfaceYZ,...
    grid.cell_surfaceXZ,grid.cell_surfaceYZ,...
    grid.cell_surfaceXY,grid.cell_surfaceXY};

grid.cell(i).faces.space_step = {grid.space_step_list(2),grid.space_step_list(1),grid.space_step_list(2),grid.space_step_list(1),grid.space_step_list(3),grid.space_step_list(3)};

%end
end


case 'darcy'
    disp('darcy law')

% add 
% 0. viscosity
% 1. K tensor of each cell (xx, yy, zz)
% 1. surface areas of each cell face
% 2. distance between cell centroid and each face centre

for i=1:length(grid.cell)

for j=1:length(grid.cell(i).faces)

% viscosity
grid.cell(i).mu = {3e-3,3e-3,3e-3,3e-3,3e-3,3e-3};  % Pa.s

% cell K in main directions
grid.cell(i).faces.K = {K(2,2),K(1,1),K(2,2),K(1,1),K(3,3),K(3,3)};
    
% cell surface areas
grid.cell(i).faces.A = {grid.cell_surfaceXZ,grid.cell_surfaceYZ,...
    grid.cell_surfaceXZ,grid.cell_surfaceYZ,...
    grid.cell_surfaceXY,grid.cell_surfaceXY};

% distance cell centroid and face centre
grid.cell(i).faces.space_step = {grid.space_step_list(2)/2,grid.space_step_list(1)/2,grid.space_step_list(2)/2,grid.space_step_list(1)/2,grid.space_step_list(3)/2,grid.space_step_list(3)/2};

%grid.cell(i).faces.T = {0,0,0,0,0,0}; %transmissibilities to be filled later


end
end
% orientation: (south, east, north, west, bottom, top)

end

%% Set up boundary struct and specific boundary conditions

% orientation: (south, east, north, west, bottom, top)
BC_struct = createBC_struct();

N_s=1; % SOUTH
N_e=2; % EAST
N_n=3; % NORTH
N_w=4; % WEST
N_b=5; % BOTTOM
N_t=6; % TOP
N_bi = 7; % bottom - inlet
N_bo = 8; % bottom - outlet


%% Impose specific conditions

% if we want periodic BC, replace "P" on respective boundary 

% Neumann (all - general)
a = 1; b = 0; c = 0;
BC_struct{1,N_w}{1,1} = a; BC_struct{1,N_w}{2,1} = b; BC_struct{1,N_w}{3,1} = c; % south
BC_struct{1,N_e}{1,1} = a; BC_struct{1,N_e}{2,1} = b; BC_struct{1,N_e}{3,1} = c; % north
BC_struct{1,N_s}{1,1} = a; BC_struct{1,N_s}{2,1} = b; BC_struct{1,N_s}{3,1} = c; % bottom
BC_struct{1,N_n}{1,1} = a; BC_struct{1,N_n}{2,1} = b; BC_struct{1,N_n}{3,1} = c; % top

BC_struct{1,N_b}{1,1} = a; BC_struct{1,N_b}{2,1} = b; BC_struct{1,N_b}{3,1} = c; 
BC_struct{1,N_t}{1,1} = a; BC_struct{1,N_t}{2,1} = b; BC_struct{1,N_t}{3,1} = c; 

% Dirichlet inlet (spiral artery entrance - bottom cells)
a = 0; b = 1; c = inlet_P; % Pin = 2506 Pa (term placenta); = 490 Pa (20weeks)
BC_struct{1,N_bi}{1,1} = a; BC_struct{1,N_bi}{2,1} = b; BC_struct{1,N_bi}{3,1} = c; 

% Dirichlet outlet (decidual veins - bottom cells)
a = 0; b = 1; c = outlet_P;
BC_struct{1,N_bo}{1,1} = a; BC_struct{1,N_bo}{2,1} = b; BC_struct{1,N_bo}{3,1} = c; 

% Dirichlet inlet (bottom)
% a = 0; b = 1; c = 1;
% BC_struct{1,N_b}{1,1} = a; BC_struct{1,N_b}{2,1} = b; BC_struct{1,N_b}{3,1} = c; 

% Dirichlet outlet (top)
% a = 0; b = 1; c = 0;
% BC_struct{1,N_t}{1,1} = a; BC_struct{1,N_t}{2,1} = b; BC_struct{1,N_t}{3,1} = c; 

% Neumann (south, north, west and east)
% a = 1; b = 0; c = 0;
% BC_struct{1,N_w}{1,1} = a; BC_struct{1,N_w}{2,1} = b; BC_struct{1,N_w}{3,1} = c; % south
% BC_struct{1,N_e}{1,1} = a; BC_struct{1,N_e}{2,1} = b; BC_struct{1,N_e}{3,1} = c; % north
% BC_struct{1,N_s}{1,1} = a; BC_struct{1,N_s}{2,1} = b; BC_struct{1,N_s}{3,1} = c; % bottom
% BC_struct{1,N_n}{1,1} = a; BC_struct{1,N_n}{2,1} = b; BC_struct{1,N_n}{3,1} = c; % top

% Periodic
%BC_struct{1,N_b} = 'P';
%BC_struct{1,N_t} = 'P';
%BC_struct{1,N_s} = 'P';
%BC_struct{1,N_n} = 'P';


%% Define cells with inlet/outlet conditions (placentome)
min_z = min(grid.cell_centroid(:,3));
min_bz = min_z-grid.space_step_list(3)/2;

s_x = 2; s_y = 2; % size of each boundary
c_x = grid.space_step_list(1); c_y = grid.space_step_list(2); % size of each cell 

nb_cellx = s_x/c_x; nb_celly = s_y/c_y; nb_cellt = nb_cellx*nb_celly;


% location of inlet
li = [network.center_position(1) network.center_position(2) min_bz];

corner_position=network.center_position-tissue_dimensions/2; 

lo1 = [corner_position(1)+tissue_dimensions(1)*0.15...
    network.center_position(2) min_bz];
lo2 = [corner_position(1)+tissue_dimensions(1)-tissue_dimensions(1)*0.15...
    network.center_position(2) min_bz];

% find cells closer to locations which fullfil desired size
cell_pos = find(grid.cell_centroid(:,3) == min_z);
cells = grid.cell_centroid(cell_pos,:);

d = sqrt((li(1)-cells(:,1)).^2 + (li(2)-cells(:,2)).^2 + (li(3)-cells(:,3)).^2);
[val,i_li] = mink(d,floor(nb_cellt));
d = sqrt((lo1(1)-cells(:,1)).^2 + (lo1(2)-cells(:,2)).^2 + (lo1(3)-cells(:,3)).^2);
[val,i_lo1] = mink(d,floor(nb_cellt));
d = sqrt((lo2(1)-cells(:,1)).^2 + (lo2(2)-cells(:,2)).^2 + (lo2(3)-cells(:,3)).^2);
[val,i_lo2] = mink(d,floor(nb_cellt));


%%
 figure; hold on; 
 scatter3(li(1),li(2),li(3)); scatter3(lo1(1),lo1(2),lo1(3)); scatter3(lo2(1),lo2(2),lo2(3));
 gpatch(meshStruct.faces, meshStruct.nodes); 
xlabel('x')
ylabel('y')
zlabel('z')


%% Assign BC struct values with faces
% orientation: (south, east, north, west, bottom, top)
for i=1:length(grid.cell)

% create initial BC structure for all faces
grid.cell(i).faces.BC = {{0;0;0},{0;0;0},{0;0;0},{0;0;0},{0;0;0},{0;0;0}};  

neigh = grid.cell(i).neigh;

for k=1:length(neigh)

    if ismember(i,i_li) == 1 && neigh(k) == -1 && k == 5 % if it is a boundary cell at bottom (inlet)

grid.cell(i).faces.BC{1,k} = BC_struct{1,7}; 

    elseif ismember(i,i_lo1) == 1 && neigh(k) == -1 && k == 5 % if it is a boundary cell at bottom (outlet1)

grid.cell(i).faces.BC{1,k} = BC_struct{1,8}; 

    elseif ismember(i,i_lo2) == 1 && neigh(k) == -1 && k == 5 % if it is a boundary cell at bottom (outlet2)

grid.cell(i).faces.BC{1,k} = BC_struct{1,8}; 

    elseif neigh(k) == -1

grid.cell(i).faces.BC{1,k} = BC_struct{1,k}; 

    end
end
end


%% Add ghost cells to boundary cells

actual_length = length(grid.cell);
grid_old = grid;

for i=1:length(grid.cell)

neigh = grid.cell(i).neigh;

for k=1:length(neigh)

    if neigh(k) == -1

grid.cell(length(grid.cell)+1).neigh = i;
grid.cell(i).neigh(k) = length(grid.cell);

    end

end
end


%% Calculate matrices: diffusion matrix, BC matrix and RHS
[M_bc, RHS,M_diff,tran] = BC_3D(law,grid,nb_cell_length,nb_cell_height,nb_cell_width);


%% SOLVE stationary or time-dependent system

%% 

switch solving

case 'stationary'
Mt = M_diff + M_bc; % matrix of coefficients for the PDE

c = Mt\RHS;
c = c(1:actual_length);


case 'time-dependent'

% time-dependent variables
dt = 0.1; % time step
final_t = 50;
alfa = 1; %coefficient of transient term

% define initial conditions
c_old = zeros(length(grid.cell),1);

% solve system
for t=dt:dt:final_t
   [M_trans, RHS_trans] = transient_struct(grid,c_old,dt,alfa);
    Mt = M_trans-M_diff+M_bc;
    RHSt = RHS_trans+RHS;
    c = Mt\RHSt;    
    c_old = c;
end


%c = c(1:actual_length);

end

%% compute velocities
% u_ij = -t_ij*(p_j - p_i)
q = zeros(length(c),6);
u = zeros(length(c),6);
u_cell = zeros(length(c),3); 

for i=1:length(grid_old.cell)
     neigh = grid_old.cell(i).neigh;

for k=1:length(neigh)
    if neigh(k) ~= -1 % for actual cells, not ghost cells
q(i,k) = -tran(i,k)*(c(neigh(k))-c(i)); % flux discharge in each direction
u(i,k) = 0.5*q(i,k)*(1/cell2mat(grid.cell(i).faces.A(k))); % contribution to cell velocity
    end
 end
 end

% compute average velocity in each cell in x,y,z directions
%(south, east, north, west, bottom, top)
for i=1:length(u)
u_cell(i,2) = -u(i,1)+u(i,3);
u_cell(i,1) = u(i,2)-u(i,4);
u_cell(i,3) = -u(i,5)+u(i,6);
end

% Average velocity magnitude in the cell
u_mag = sqrt(u_cell(:,1).*u_cell(:,1) + u_cell(:,2).*u_cell(:,2) + u_cell(:,3).*u_cell(:,3));

% normalize
for i=1:length(u_cell)
un_cell(i,:) = u_cell(i,:)/norm(u_cell(i,:));
end




%% plot pressure
p = reshape(c,nb_cell_length,nb_cell_height,nb_cell_width);
p = full(p);

 figure
  h=slice(p,nb_cell_length/2,nb_cell_height/2,nb_cell_width/2); hold on;
  set(h, 'EdgeColor','none', 'FaceColor','interp');
 a = colorbar; clim([100 180]);
 %plot_tree(tree)
 a.Label.String = 'Pressure (Pa)';
 xlabel('x')
 ylabel('y')
 zlabel('z')

%% plot velocity magnitude
u_mag2 = reshape(u_mag,nb_cell_length,nb_cell_height,nb_cell_width);
u_mag2 = full(u_mag2);
mat = [grid.cell_centroid(:,1) grid.cell_centroid(:,2) grid.cell_centroid(:,3)...
    u_cell(:,1) u_cell(:,2) u_cell(:,3)];

 figure; hold on;
  h=slice(u_mag2,nb_cell_length/2,nb_cell_height/2,nb_cell_width/2);
  set(h, 'EdgeColor','none', 'FaceColor','interp');
  %quiver3(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),4)

 a = colorbar; %clim([0 200]);
 a.Label.String = 'Velocity magnitude (mm^2/s)';
 xlabel('x')
 ylabel('y')
 zlabel('z')


%% plot velocity vectors

mat = [grid.cell_centroid(:,1) grid.cell_centroid(:,2) grid.cell_centroid(:,3)...
    u_cell(:,1) u_cell(:,2) u_cell(:,3)];

figure; hold on;
quiver3(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),6,'r')
axis equal

%% PLOT -- mute this bit if no Gibbon installed

var = u_mag;

% % color codes
  numpoints = ceil(length(var));
  cmap = parula(numpoints);
  range_c = linspace(min(var), max(var),numpoints)';
% % 
  CF = [];
% % 
  for i=1:length(var)
% % 
  [m,I] = min(abs(range_c(:)-var(i)));  
  CF = [CF; cmap(I,:)];
  end
% % 
  meshStruct.elementData = var;
% % 
  defaultOptionStruct.numSLiceSteps=25; %Number of slice steps
  defaultOptionStruct.cMap=CF; %Colormap used (if empty it is based of number of element material types
  defaultOptionStruct.faceAlpha1=0.2; %Alpha level for boundary surface
  defaultOptionStruct.faceAlpha2=1; %Alpha level for mesh elements
  defaultOptionStruct.lightWeightPlot=1; %Option to only plot element outer boundaries to create a more lightweigth plot
  defaultOptionStruct.cMap=CF;
% % 

figure; hold on;
  meshView(meshStruct,defaultOptionStruct); plot_tree(tree)
%quiver3(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),6,'r')

% % 
 % axisGeom(gca,15);
  data_range = 50;
  colormap(parula(data_range));
  colorbar();
  a = colorbar; clim([0 50]);
  a.Label.String = 'Velocity (mm2/s)';
% % 
% % hold on; plot_tree(tree)
  drawnow;
% 
