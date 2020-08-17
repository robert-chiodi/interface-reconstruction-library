% This file is part of the Interface Reconstruction Library (IRL),
% a library for interface reconstruction and computational geometry operations.
%
% Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at https://mozilla.org/MPL/2.0/.

clear all; clf; close all; 

% Read in the mesh data
viz_directory = '/path/to/viz/folder';
meshfileID = fopen(strcat(viz_directory,'/mesh'),'r');
nx = fscanf(meshfileID,'%i',1);
x = fscanf(meshfileID,'%E',nx+1);
ny = fscanf(meshfileID,'%i',1);
y = fscanf(meshfileID,'%E',ny+1);
nz = fscanf(meshfileID,'%i',1);
z = fscanf(meshfileID,'%E',nz+1);
fclose(meshfileID);

xm(1:nx) = 0.0;
ym(1:ny) = 0.0;
zm(1:nz) = 0.0;
for i = 1:nx
   xm(i) = 0.5*(x(i)+x(i+1)); 
end
for j = 1:ny
   ym(j) = 0.5*(y(j)+y(j+1)); 
end
for k = 1:nz
   zm(k) = 0.5*(z(k)+z(k+1)); 
end

directory = dir([viz_directory, '/vizfile_*']);
number_of_viz_files = size(directory,1);
liquid_volume_fraction(1:nx,1:ny) = 0.0;
for iteration = 0:number_of_viz_files-1
    filename = strcat(viz_directory,strcat('vizfile_',int2str(iteration)));
    datafileID = fopen(filename);
    liquid_volume_fraction(1:nx,1:ny) = reshape(fscanf(datafileID,'%E',nx*ny*nz),[ny,nx]);
    fclose(datafileID);
figure(1)
    pcolor(xm,ym,liquid_volume_fraction);
    colorbar;
    caxis([0.0 1.0]);
    axis([x(1) x(nx+1) y(1) y(ny+1)]);
   formatSpec = strcat(viz_directory,'/VOF_%d');
   pname = sprintf(formatSpec,iteration);
   print(pname,'-dpng')
end
