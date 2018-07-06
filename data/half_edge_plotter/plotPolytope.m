%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function [] = plotPolytope( nfaces, nvert_on_faces, vertices_on_faces, ref_point, figure_number, color)

figure(figure_number)
hold on;
for n = 1:nfaces
   plotFace(nvert_on_faces(n),vertices_on_faces(:,:,n), ref_point, figure_number, color); 
end

end

