
function [M_trans, RHS_trans] = transient_struct(grid,c_old, dt, alfa)

% create BC matrix and RHS structures
M_trans = sparse(length(grid.cell),length(grid.cell));
RHS_trans = zeros(length(grid.cell),1);

%% UPDATE TRANSIENT TERMS
for i=1:length(grid.cell)

if length(grid.cell(i).neigh) > 1  % internal cells

M_trans(i,i) = alfa/dt;

% add value to RHS (row of internal cells)
RHS_trans(i,1) =  alfa.*c_old(i)/dt;

end  
end



end