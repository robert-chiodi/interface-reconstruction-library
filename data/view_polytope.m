%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

vertex(1:3,1) = [ -1.5, -1.5, -0.5 ];
neighbor(1:3,1) = [5,2,1];
vertex(1:3,2) = [ -1.5, 1.5, -0.5 ];
neighbor(1:3,2) = [0,2,3];
vertex(1:3,3) = [ -1.5, 1.5, 0.5 ];
neighbor(1:3,3) = [4,1,0];
vertex(1:3,4) = [ -0.5, 1.5, -0.5 ];
neighbor(1:3,4) = [1,4,5];
vertex(1:3,5) = [ -0.5, 1.5, 0.166666666666667 ];
neighbor(1:3,5) = [2,5,3];
vertex(1:3,6) = [ -0.5, -0.5, -0.5 ];
neighbor(1:3,6) = [0,3,4];

neighbor = neighbor + 1;

edges_visited(1:3,1:length(vertex(1,:))) = 0;
edge_index(1:2) = 0;
figure(1)
for v = 1:length(vertex(1,:))
    for n = 1:3

        if(edges_visited(n,v) == 0)
            edge_index(1) = v;
            edge_index(2) = neighbor(n,v);
            for neigh = 1:3
               if(neighbor(neigh,edge_index(2)) == edge_index(1))
                  break; 
               end
            end
            edge_index(1) = edge_index(2);
            next_index  = mod((neigh),3)+1;
            edge_index(2) = neighbor(next_index,edge_index(1));
            edges_visited(next_index, edge_index(1)) = 1;
            
            while(edge_index(1) ~= v)
            tri_vertex(:,1) = vertex(:,v);
            tri_vertex(:,2) = vertex(:,edge_index(1));
            tri_vertex(:,3) = vertex(:,edge_index(2));           
            
            fill3(tri_vertex(1,:),tri_vertex(2,:),tri_vertex(3,:),'r');
            hold on;
            
            for neigh = 1:3
               if(neighbor(neigh,edge_index(2)) == edge_index(1))
                  break; 
               end
            end
            edge_index(1) = edge_index(2);
            next_index  = mod((neigh),3)+1;
            edge_index(2) = neighbor(next_index,edge_index(1));
            edges_visited(next_index, edge_index(1)) = 1;
                
            end

        end
    end
    
end
hold off;
xlabel('x');
ylabel('y');
zlabel('z');
