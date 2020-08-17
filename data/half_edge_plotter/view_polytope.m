%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

%clear all;
%close all;
% Load in the polytope half-edges
run('polytope.m');

graphics_toolkit("gnuplot")

ref_point = [-1.0 -0.5 0.5];

% Plot an individual face
face = 0;
% for face = 1:1
%     face_index_to_plot = face;
%     plotFace(nvert_for_face(face_index_to_plot), ...
%         vert_on_face(:,:,face_index_to_plot), ref_point, face);
% end
 
% Plot the entire polytope
plotPolytope(nfaces, nvert_for_face, vert_on_face, ref_point, face+3,"g");

run('polytope2.m');
plotPolytope(nfaces, nvert_for_face, vert_on_face, ref_point, face+4,"r");
