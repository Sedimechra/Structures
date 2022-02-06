angle = 0;
rotation = [cosd(angle) -sind(angle) 0 0; sind(angle) cosd(angle) 0 0; 0 0 cosd(angle) -sind(angle); 0 0 sind(angle) cosd(angle)];
forces2 = [7.48; 0; 0; -200];
g = transpose(rotation) * forces2