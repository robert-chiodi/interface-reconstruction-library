%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

% This will create the lookup table for cutting a PLIC polygon surface
% by a plane. The lookup table will be used to truncate multiple planes
% in a cell to just have the final geometry

% This lookup table is only valid for starting polygons
% that start with 3-6 vertices and can end with 3-7 vertices.

clear all;

d(1:7) = 0.0;
intersection_points(1:64) = int8(-1); % Number of intersection pts between plane and polygon
npts_in_new_poly(1:64,1:6) = int8(-1); % Number of points in polygon after cutting
edge_vertex_1(1:64,1:2) = int8(0); % Vertex 1 on intersection edges (to interpolate for point)
edge_vertex_2(1:64,1:2) = int8(0); % Vertex 2 on intersection edges (to interpolate for point)
verts_in_new_poly(1:64,1:7) = int8(0); % List of vertices in new polygon (vertices 8,9 are intersection points)

write_out_c_lookup = true;

% These are the signs of the 6 vertices for
% our max 6 vertex polygon.
for v6 = -1:2:1
    for v5 = -1:2:1
        for v4 = -1:2:1
            for v3 = -1:2:1
                for v2 = -1:2:1
                    for v1 = -1:2:1                       
                        np = 6;
                        d(1) = v1;
                        d(2) = v2;
                        d(3) = v3;
                        d(4) = v4;
                        d(5) = v5;
                        d(6) = v6;
                        d(np+1:7) = d(1);
                        ca = 1+int8(0.5+0.51*sign(d(1)))+...
                             2*int8(0.5+0.51*sign(d(2)))+...
                             4*int8(0.5+0.51*sign(d(3)))+...
                             8*int8(0.5+0.51*sign(d(4)))+...
                             16*int8(0.5+0.51*sign(d(5)))+...
                             32*int8(0.5+0.51*sign(d(6))); 
                        % If has >2 sign changes, cycle
                        % because this cant happen when cut by a plane
                        sumchange = 0;
                        for i = 1:np
                            if(d(i)*d(i+1) < 0) 
                                sumchange = sumchange + 1;
                            end
                        end
                        if(ca == 61)
                           fprintf('Number of changes: %d \n',sumchange);     
                        end
                        if(sumchange > 2)
                            continue;
                        end
                        
                         % Find crossing and record the vertices involved
                             npts = 0;
                             for i = 1:np
                                 if(d(i)*d(i+1)<0)
                                    npts = npts + 1;
                                    edge_vertex_1(ca,npts) = i;
                                    edge_vertex_2(ca,npts) = i+1;
                                    % The np+1 vertex is actually vertex 1
                                    if(edge_vertex_2(ca,npts) == np+1)
                                        edge_vertex_2(ca,npts) = 1;
                                    end
                                 end
                             end
                             intersection_points(ca) = npts;

                             % Now using knowledge of crossings,
                             % record vertices involved
                             % in the new polygon
                             % (which are all the negative ones
                             % and the crossings)
                             cnt = 0;
                             v = 1;
                             for i = 1:np
                                 % Crossing but on positive side, need the
                                 % intersection point
                                 if(d(i)*d(i+1)<0 && d(i)>0)                          
                                     verts_in_new_poly(ca,v) = 8+cnt; 
                                     cnt = cnt + 1;
                                     v = v + 1;             
                                 end
                                 % Negative, so need this point
                                 if(d(i) < 0) 
                                    verts_in_new_poly(ca,v) = i; 
                                    v = v + 1;
                                 end
                                 % Crossing but on negative side
                                 if(d(i)*d(i+1)<0 && d(i)<0)
                                    verts_in_new_poly(ca,v) = 8+cnt; 
                                    v = v + 1;
                                    cnt = cnt + 1;
                                 end
                             end
                    end
                end
            end
        end
    end
end

for v6 = -1:2:1
    for v5 = -1:2:1
        for v4 = -1:2:1
            for v3 = -1:2:1
                for v2 = -1:2:1
                    for v1 = -1:2:1     
                        for np = 3:6;
                            d(1) = v1;
                            d(2) = v2;
                            d(3) = v3;
                            d(4) = v4;
                            d(5) = v5;
                            d(6) = v6;
                            d(np+1:7) = d(1);    
                            ca = 1+int8(0.5+0.5*sign(d(1)))+...
                             2*int8(0.5+0.5*sign(d(2)))+...
                             4*int8(0.5+0.5*sign(d(3)))+...
                             8*int8(0.5+0.5*sign(d(4)))+...
                             16*int8(0.5+0.5*sign(d(5)))+...
                             32*int8(0.5+0.5*sign(d(6))); 
                        % If has >2 sign changes, cycle
                        % because this cant happen when cut by a plane
                        sumchange = 0;
                        positive_in_range = 0;
                        for i = 1:np
                            if(d(i)*d(i+1) < 0) 
                                sumchange = sumchange + 1;
                            end
                            if(d(i) > 0)
                                positive_in_range = positive_in_range+1;
                            end
                        end
                        if(sumchange > 2)
                            continue;
                        end

                            
                            npts_in_new_poly(ca,np) = np - ...
                            positive_in_range + sumchange;
                            
                        
                        end
                    end
                end
            end
        end
    end
end
                        
                        
                        
                        


%  i = 0;
%  fprintf('\n integer, parameter, dimension(64) :: pc_npts =(/&\n')
%  for ca2 = 1:64
%      i = i + 1;
%              fprintf('%d,',pc_npts(ca2,6));
%      if(i >= 32)
%         fprintf('&\n')
%         i = 0;
%      end
%  end
%  fprintf('/)')
%  
%   i = 0;
%   fprintf('\n integer, parameter, dimension(6,64) :: pc_v1 =reshape((/&\n')
%  for ca2 = 1:64
%              for v2 = 1:6
%                       i = i + 1;
%                 fprintf('%d,',pc_v1(ca2,6,v2));    
%              end
%      if(i >= 32)
%         fprintf('&\n')
%         i = 0;
%      end
%  end
%  fprintf('/),(/6,64/))')
%  
%   i = 0;
%    fprintf('\n integer, parameter, dimension(6,64) :: pc_v2 =reshape((/&\n')
%  for ca2 = 1:64
%              for v2 = 1:6
%                       i = i + 1;
%                 fprintf('%d,',pc_v2(ca2,6,v2));    
%              end
%      if(i >= 32)
%         fprintf('&\n')
%         i = 0;
%      end
%  end
%  fprintf('/),(/6,64/))')
%  
%   i = 0;
%     fprintf('\n integer, parameter, dimension(9,64) :: pc_vert =reshape((/&\n')
%  for ca2 = 1:64
%              for v2 = 1:9
%                       i = i + 1;
%                 fprintf('%d,',pc_vert(ca2,6,v2));    
%              end
%      if(i >= 32)
%         fprintf('&\n')
%         i = 0;
%      end
%  end
%  fprintf('/),(/9,64/))')
%  
%   i = 0;
%  fprintf('\n integer, parameter, dimension(6,64) :: pc_newpts =reshape((/&\n')
%  for ca2 = 1:64
%          for v = 1:6
%                   i = i + 1;
%              fprintf('%d,',pc_newpts(ca2,v));
%          end
%      if(i >= 32)
%         fprintf('&\n')
%         i = 0;
%      end
%  end
%  fprintf('/),(/6,64/))')
    

if(write_out_c_lookup)
    % Print functions for C lookup table
    i = 0;
     fprintf('\n constexpr int number_of_plane_intersections_with_polygon[64] = { \n')
     for ca2 = 1:64
         i = i + 1;
                 fprintf('%d, ',intersection_points(ca2));
         if(i >= 16)
            fprintf('\n')
            i = 0;
         end
     end
     fprintf('};')

     i = 0;
     fprintf('\n constexpr int number_of_points_after_intersecting_plane[64][6] = { \n')
     for ca2 = 1:64
         fprintf('{');
         for v=1:6
            if(v<6)
                fprintf('%d, ',npts_in_new_poly(ca2,v));
            end
            if(v==6)
                fprintf('%d ',npts_in_new_poly(ca2,v));
            end
         end
         fprintf('},\n');
     end
     fprintf('};\n')

     fprintf('\n constexpr int plane_cut_vertices[64][2][2] = { \n')
     for ca2 = 1:64
         fprintf('{');
         % Vertex v1
             fprintf('{');
            for v = 1:2
                if(v < 2)
                    fprintf('%d, ',edge_vertex_1(ca2,v)-1);
                else
                    fprintf('%d ',edge_vertex_1(ca2,v)-1);
                end
            end
         fprintf('},\n'); 
         fprintf('{');
         % Vertex v2
            for v = 1:2
                if(v < 2)
                    fprintf('%d, ',edge_vertex_2(ca2,v)-1);
                else
                    fprintf('%d ',edge_vertex_2(ca2,v)-1);
                end
            end
         fprintf('}');
         fprintf('}, \n');
     end
     fprintf('};\n')

     i = 0;
     fprintf('\n constexpr int vertex_to_add[64][7] = { \n')
     for ca2 = 1:64
         fprintf('{');
            for v=1:7
                if(v<7)
                    fprintf('%d, ',verts_in_new_poly(ca2,v)-1);
                end
                if(v==7)
                    fprintf('%d ',verts_in_new_poly(ca2,v)-1);
                end
            end
         fprintf('},\n');
     end
     fprintf('}; \n')
end


     
