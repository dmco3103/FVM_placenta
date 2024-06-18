function grid=generateGrid2(grid_center_position,dimension_list, nb_cell_list)

%% Prep for grid setup

%%
corner_position=grid_center_position-dimension_list/2; %[width/2 height/2];
space_step_list=dimension_list./nb_cell_list;
% dX=width/nb_cell_width;
% dY=height/nb_cell_height;
cell_position_list=[];

for k=1:nb_cell_list(3) % z
    for j=1:nb_cell_list(2) % y 
        for i=1:nb_cell_list(1) % x
            l=i+(j-1)*nb_cell_list(1)+(k-1)*nb_cell_list(1)*nb_cell_list(2);
            position=[];
            indices=[i j k];
            for h=1:length(dimension_list)
                position=[position corner_position(h)+(1/2+(indices(h)-1))*space_step_list(h)];
            end

% DEFINE, each face as internal -I- or boundary -B-  AND define neighbours of each face       
%neigh = [s e n w b t];
if i == 1
    w = -1; type_w = 'B';
    e = l+1; type_e = 'I';
elseif i == nb_cell_list(1)
    e = -1; type_e = 'B';
    w = l-1; type_w = 'I';
else
    e = l+1; w = l-1; type_e = 'I'; type_w = 'I';
end

if j == 1
    s = -1; type_s = 'B';
    n = l+nb_cell_list(2); type_n = 'I';
elseif j == nb_cell_list(2)
    n = -1; type_n = 'B';
    s = l-nb_cell_list(2); type_s = 'I';
    else
    s = l-nb_cell_list(2); n = l+nb_cell_list(2); type_s = 'I'; type_n = 'I';
end

if k == 1
    b = -1; type_b = 'B';
    t = l + nb_cell_list(1)*nb_cell_list(2); type_t = 'I';
elseif k == nb_cell_list(3)
    t = -1; type_t = 'B';
    b = l - nb_cell_list(1)*nb_cell_list(2); type_b = 'I';
else
    t = l + nb_cell_list(1)*nb_cell_list(2); type_t = 'I';
    b = l - nb_cell_list(1)*nb_cell_list(2); type_b = 'I';
end


            tissue.cell(l).centroid=position;
            tissue.cell(l).neigh = [s e n w b t];
            tissue.cell(l).faces.name = {'s','e','n','w','b','t'};
            tissue.cell(l).faces.type = {type_s,type_e,type_n,type_w,type_b,type_t};

            cell_position_list=[cell_position_list; position ];
        end
    end
end

tissue.cell_centroid=cell_position_list;
tissue.dimension_list=dimension_list;
tissue.nb_cell_list=nb_cell_list;
tissue.space_step_list=space_step_list;

% cell areas & volume
tissue.cell_surfaceXY = tissue.space_step_list(1)*tissue.space_step_list(2);
tissue.cell_surfaceXZ = tissue.space_step_list(1)*tissue.space_step_list(3);
tissue.cell_surfaceYZ = tissue.space_step_list(2)*tissue.space_step_list(3);

tissue.cell_volume=space_step_list(1)*space_step_list(2)*space_step_list(3);

grid = tissue;

end

