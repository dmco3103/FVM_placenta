function BC = createBC_struct()

% Creates a boundary condition structure (2D for now): (S E N W) 

% structure: a*(grad(T)) + b*T = c
% Dirichlet: a=0;b~=0
% Neumann: a~=0; b=0
% Robin: a~=0; b~=0

% STRUCTURE: list of cells with 3x1 dimensions for each position

% BC = [S  E  N  W]
%      [a  a  a  a]
%      [b  b  b  b]
%      [c  c  c  c]



% Define the top, bottom, right, and left boundary conditions
% All default with a = 1 (Neumann)

BC = {{1;0;0},{1;0;0},{1;0;0},{1;0;0}};  


end