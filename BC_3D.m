
function [M_bc, RHS,M_diff,tran] = BC_3D(law,grid,nb_cell_length,nb_cell_height,nb_cell_width)

M_diff = sparse(length(grid.cell),length(grid.cell));

% create BC matrix and RHS
M_bc = sparse(length(grid.cell),length(grid.cell));
RHS = zeros(length(grid.cell),1);

tran = [];

%% ALL CONDITIONS
for i=1:length(grid.cell)

if length(grid.cell(i).neigh) > 1 % actual cells (no neighbours)

% get neighbours (south, east, north, west, bottom, top) - will correspond to column in matrix M
neigh = grid.cell(i).neigh;


switch law

case 'diffusion'

for k=1:length(neigh)

% DIFFUSION contribution from internal cells to their own columns and diagonal entry
        M_diff(i,neigh(k)) = cell2mat(grid.cell(i).faces.D(k))* ...
            cell2mat(grid.cell(i).faces.A(k)) / ...
            cell2mat(grid.cell(i).faces.space_step(k)); 

        M_diff(i,i) = M_diff(i,i) - M_diff(i,neigh(k));
end

case 'darcy'

for k=1:length(neigh)

    % half-t (main cell)
   % t = (K_a * S) / mu*r_a ; 
t1 = (cell2mat(grid.cell(i).faces.K(k)) * cell2mat(grid.cell(i).faces.A(k))) / ...
    (cell2mat(grid.cell(i).mu(k)) * cell2mat(grid.cell(i).faces.space_step(k)));

if neigh(k) > size(grid.cell_centroid,1) % if neighbour is ghost cell
    % half-t (neighbour)
t2 = (cell2mat(grid.cell(i).faces.K(k)) * cell2mat(grid.cell(i).faces.A(k))) / ...
    (cell2mat(grid.cell(i).mu(k)) * cell2mat(grid.cell(i).faces.space_step(k)));

else
%south, east, north, west, bottom, top
if k == 1
    c_P = 3;
elseif k == 2
    c_P = 4;
elseif k == 3
    c_P = 1;
elseif k == 4
    c_P = 2;
elseif k == 5
    c_P = 6;
elseif k == 6
    c_P = 5;
end

t2 = (cell2mat(grid.cell(neigh(k)).faces.K(c_P)) * cell2mat(grid.cell(i).faces.A(k))) / ...
    (cell2mat(grid.cell(neigh(k)).mu(c_P)) * cell2mat(grid.cell(i).faces.space_step(k)));
end
    
    % compute total t and store
T = (1/t1 + 1/t2)^-1;    
tran(i,k) = T;

% Flux contribution from internal cells to their own columns and diagonal entry
        M_diff(i,neigh(k)) = -T; 

        M_diff(i,i) = M_diff(i,i) - M_diff(i,neigh(k));
end 

%grid.cell(i).faces.T = {tran(1),tran(2),tran(3),tran(4),tran(5),tran(6)};


end


% get face type - interior or boundary
t = grid.cell(i).faces.type;

for k=1:length(t)

  %  NON-PERIODIC CONDITIONS
    if ismember(t(k),'B') == 1 && length(grid.cell(i).faces.BC{1,k}) > 1  % if it is a boundary, we apply BC fluxes to this cell and its ghost //  % for non-periodic conditions

% if not periodic, get BC (a,b,c) values
bc = cell2mat(grid.cell(i).faces.BC{1,k});

orient = cell2mat(grid.cell(i).faces.name(k));

switch orient
    case 's'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
     -(bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
     -(bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
RHS(grid.cell(i).neigh(k),1) =  -bc(3);


    case 'e'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
    (bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
    (bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
RHS(grid.cell(i).neigh(k),1) = bc(3);


    case 'n'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
    (bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
    (bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
RHS(grid.cell(i).neigh(k),1) = bc(3);
       

    case 'w'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
     -(bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
     -(bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
RHS(grid.cell(i).neigh(k),1) =  -bc(3);


 case 'b'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
     -(bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
     -(bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
RHS(grid.cell(i).neigh(k),1) =  -bc(3);


case 't'

% add coefficient to ghost cell itself and influence on actual cell
M_bc(grid.cell(i).neigh(k),grid.cell(i).neigh(k)) = ...
    (bc(2)/2 + bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

M_bc(grid.cell(i).neigh(k),i) = ...
    (bc(2)/2 - bc(1)/cell2mat(grid.cell(i).faces.space_step(k)));

% add value to RHS
    RHS(grid.cell(i).neigh(k),1) = bc(3);

end

% PERIODIC CONDITION
    elseif ismember(t(k),'B') == 1 && length(grid.cell(i).faces.BC{1,k}) == 1  % if it is a boundary, we apply BC fluxes to this cell and its ghost //  % for non-periodic conditions
% 
 orient = cell2mat(grid.cell(i).faces.name(k));
% 
 switch orient
% 
     case 's'
 neighbour = grid.cell(i).neigh(k);
 opposite = i+nb_cell_length*(nb_cell_height-1);
 opposite_neighbour = grid.cell(opposite).neigh(3);
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = 1;
 M_bc(neighbour,opposite_neighbour) = -1;
% 
% % add value to RHS
 RHS(neighbour,1) =  0;
% 
% 
     case 'e'
% 
% % add coefficient to ghost cell itself and influence on actual cell
 neighbour = grid.cell(i).neigh(k);
 opposite = i - (nb_cell_length-1);
 opposite_neighbour = grid.cell(opposite).neigh(4);
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = -1;
 M_bc(neighbour,opposite_neighbour) = 1;
% 
% % add value to RHS
 RHS(neighbour,1) =  0;
% 
% 
     case 'n'
% 
% % add coefficient to ghost cell itself and influence on actual cell
 neighbour = grid.cell(i).neigh(k);
 opposite = i-nb_cell_length*(nb_cell_height-1);
 opposite_neighbour = grid.cell(opposite).neigh(1);
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = -1;
 M_bc(neighbour,opposite_neighbour) = 1;
% 
% % add value to RHS
 RHS(neighbour,1) =  0;
% 
% 
     case 'w'
% 
 neighbour = grid.cell(i).neigh(k);
 opposite = i + (nb_cell_length-1);
 opposite_neighbour = grid.cell(opposite).neigh(2);
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = 1;
 M_bc(neighbour,opposite_neighbour) = -1;
% 
% % add value to RHS
 RHS(neighbour,1) =  0;
% 
% 
% %%
  case 'b'
% 
 neighbour = grid.cell(i).neigh(k);
 opposite = i + nb_cell_length*nb_cell_height*(nb_cell_width-1);
 opposite_neighbour = grid.cell(opposite).neigh(6);
% 
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = 1;
 M_bc(neighbour,opposite_neighbour) = -1;
% 
% % add value to RHS
 RHS(neighbour,1) =  0;
% 
% 
 case 't'
% 
% % add coefficient to ghost cell itself and influence on actual cell
 neighbour = grid.cell(i).neigh(k);
 opposite = i - nb_cell_length*nb_cell_height*(nb_cell_width-1);
 opposite_neighbour = grid.cell(opposite).neigh(5);
% 
% 
% % add coefficient to ghost cell itself and influence on actual cell
 M_bc(neighbour,neighbour) = 1;
 M_bc(neighbour,opposite) = -1;
 M_bc(neighbour,i) = -1;
 M_bc(neighbour,opposite_neighbour) = 1;
% 
% % add value to RHS
 RHS(neighbour,1) = 0;
% 
% 
% 

 end
    end
end
end
end





end


