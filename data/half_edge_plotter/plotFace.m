%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function [] = plotFace( nvert_on_face,  vertices_on_face, ref_point, figure_number, color)

%vol = 0.0; 
%for v = 3:nvert_on_face
%    for n = 1:3
%        a(n) = vertices_on_face(n,1) - ref_point(n);
%        b(n) = vertices_on_face(n,v-1) - ref_point(n);
%        c(n) = vertices_on_face(n,v) - ref_point(n);
%    end
%    vol = vol + dot(a, cross(b,c)) / 6.0;
%end

figure(figure_number)
title('Full Faces');
for n = 1:nvert_on_face-2 % Triangular for Octave gnuplot
triangle(1:3,1) = vertices_on_face(1:3,1);
triangle(1:3,2) = vertices_on_face(1:3,n+1);
triangle(1:3,3) = vertices_on_face(1:3,n+2);
h = patch(triangle(1,1:3),...
    triangle(2,1:3),...
    triangle(3,1:3),color);
 set (h, "EdgeColor", "k");    
end
hold on;
s(1:nvert_on_face) =10;
scatter3(vertices_on_face(1,1:nvert_on_face),...
    vertices_on_face(2,1:nvert_on_face),...
    vertices_on_face(3,1:nvert_on_face),s,'filled');
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1]);
pbaspect([1 1 1]);



end

