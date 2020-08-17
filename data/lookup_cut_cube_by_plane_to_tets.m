%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

clear all;
close all;

% Create unit cube
x=[0,1];
y=[0,1];
z=[0,1];
G=zeros(2,2,2);
visualize=0;
% If changing this point ordering, also need to change "faces"
% for visualization AND the edges (e1,e2) pairs when looking
% for brute force intersections
pt=zeros(4,8);
i=1; j=1; k=1; n=0;
pt(1,1) = x(2); pt(2,1) = y(1); pt(3,1) = z(1);
pt(1,2) = x(2); pt(2,2) = y(1); pt(3,2) = z(2); 
pt(1,3) = x(2); pt(2,3) = y(2); pt(3,3) = z(2);
pt(1,4) = x(2); pt(2,4) = y(2); pt(3,4) = z(1);
pt(1,5) = x(1); pt(2,5) = y(1); pt(3,5) = z(1);
pt(1,6) = x(1); pt(2,6) = y(1); pt(3,6) = z(2);
pt(1,7) = x(1); pt(2,7) = y(2); pt(3,7) = z(2);
pt(1,8) = x(1); pt(2,8) = y(2); pt(3,8) = z(1);

% Data to tabulate
table_done=zeros(1,256);
table_ninterp(1:256)=-1;
table_v(1:256,1:2,1:6)=0;

% Points to turn to vertices of tet
point_list_for_tet_building(1:14,1:3) = 0.0;
% List saying which vertex number each point in above list is
point_to_stored_vertex_map(1:14) = 0;
% Number of tets under plane in cube
number_of_tets_under_plane(1:256) = -1;
% Total number of tets, number_of_tets_under_plane + those above
total_number_of_tets(1:256) = -1;
% Vertices involved in tets
tet_vertices(1:256,1:17,1:4) = 0;

% Outer loop here - note that we have only 104 cases/256 in the presence of
% a single plane (see Boolean threshold functions)
while sum(table_done)<104

    % Create random plane
    theta1=random('unif',0,2*pi);
    theta2=random('unif',0,2*pi);
    [nplic(1),nplic(2),nplic(3)]=sph2cart(theta1,theta2,1.0);
    d=random('unif',-2,2);

    % Create 8 vertices with distance
    for n = 1:8
         pt(4,n)=nplic(1)*pt(1,n)+nplic(2)*pt(2,n)+nplic(3)*pt(3,n)-d;
        % It seems that angle calculation can cause problems - let
        % us avoid limit cases!
        if pt(4,n)>= 0
            pt(4,n)=1;
        else
            pt(4,n)=-1;
        end
    end
    
    % Hardcode a testcase
    %pt(4,1)=+8.2362171793863986E-003;
    %pt(4,2)=-1.7655694733176949E-002;
    %pt(4,3)=+2.5733971249896709E-002;
    %pt(4,4)=-1.5794066266663909E-004;
    %pt(4,5)=+8.2362171793863986E-003;
    %pt(4,6)=-1.7655694733176949E-002;
    %pt(4,7)=+2.5733971249896709E-002;
    %pt(4,8)=-1.5794066266663909E-004;
    
    % Create cube index
    iplic=1;
    for n=1:8
        if pt(4,n)>=0
            iplic=iplic+2^(n-1);
        end
    end

    % Work only if case was not treated already
    if table_done(iplic)==0
        
        % Plot out planes intersection the cube
        %visualize = 1;
        
        % Display the problem
        for n=1:8
            vertices(n,1:3)=pt(1:3,n);
        end
        faces=[1 2 3 4; 5 6 7 8; 1 2 6 5; 4 3 7 8; 1 4 8 5; 2 3 7 6];
        if visualize==1
            patch('Vertices',vertices,'Faces',faces,'EdgeColor','blue','FaceColor','none','LineWidth',1)%,'Marker','o','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',10)
            hold on
            for n=1:8
                if pt(4,n)>=0
                    scatter3(vertices(n,1),vertices(n,2),vertices(n,3),50,'red' ,'o','MarkerEdgeColor' ,'red','MarkerFaceColor','red' )
                else
                    scatter3(vertices(n,1),vertices(n,2),vertices(n,3),50,'blue','o','MarkerEdgeColor','blue','MarkerFaceColor','blue')
                end
            end
            surfx=[1.2 -0.2 -0.2 1.2]; surfy=[1.2 1.2 -0.2 -0.2]; surfz=-1/nplic(3)*(nplic(1)*surfx+nplic(2)*surfy-d);
            fill3(surfx,surfy,surfz,'red','FaceAlpha',0.4)
        end

        % We now will bruteforce create all intersections
        clear v1;
        clear v2;
        clear myp;
        ninterp=0;
        e1=1; e2=2;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=2; e2=3;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=3; e2=4;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=4; e2=1;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=5; e2=6;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=6; e2=7;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=7; e2=8;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=8; e2=5;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=1; e2=5;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=2; e2=6;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=4; e2=8;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        e1=3; e2=7;
        if pt(4,e1)*pt(4,e2)<=0
            ninterp=ninterp+1;v1(ninterp)=e1;v2(ninterp)=e2;
            mu=min(1,max(0,-pt(4,e1)/((abs(pt(4,e2)-pt(4,e1))+realmin)*sign(pt(4,e2)-pt(4,e1)))));
            myp(ninterp,1)=(1-mu)*pt(1,e1)+mu*pt(1,e2);myp(ninterp,2)=(1-mu)*pt(2,e1)+mu*pt(2,e2);myp(ninterp,3)=(1-mu)*pt(3,e1)+mu*pt(3,e2);
            if visualize==1
                scatter3(myp(ninterp,1),myp(ninterp,2),myp(ninterp,3),50,'green','d','MarkerEdgeColor','green','MarkerFaceColor','green')
            end
        end
        
        % Generate reordered polygon
        if (ninterp>0)
            % Correct order if needed
            clear myang;
            bary=1/ninterp*sum(myp,1);
            d1=myp(1,:)-bary;  d1=d1/norm(d1,2);
            d2=cross(nplic,d1);d2=d2/norm(d2,2);
            a(1)=dot(myp(1,:)-bary,d1);a(2)=dot(myp(1,:)-bary,d2);
            angref=atan2(a(2),a(1));angref=mod(angref+2*pi,2*pi);
            for n=1:ninterp
                b(1)=dot(myp(n,:)-bary,d1);b(2)=dot(myp(n,:)-bary,d2);
                myang(n)=atan2(b(2),b(1));
                myang(n)=mod(myang(n)-angref+4*pi,2*pi);
            end
            [sortedangle,myface]=sort(myang);
            % Display polygon
            if visualize==1
                patch('Vertices',myp,'Faces',myface,'FaceColor','green','EdgeColor','black','LineWidth',0.5)
            end
        end

        % Collect table
        table_done(iplic)=1;
        table_ninterp(iplic)=ninterp;
        for n=1:ninterp
            table_v(iplic,1,n)=v1(myface(n));
            table_v(iplic,2,n)=v2(myface(n));
        end
        display(sum(table_done));
        if visualize==1
            display(table_ninterp(iplic));
            display(table_v(iplic,1,:));
            display(table_v(iplic,2,:));
        end
        
        % Collect interface polygon points and negative
        % cube vertex points
        added_points = 0;
        clear point_list_for_tet_building;
        clear point_to_stored_vertex_map;
        clear avail_ptlist;
        point_list_for_tet_building(:,1:3) = 0.0;
        point_to_stored_vertex_map(:) = 0;
        avail_ptlist(1:14,1:3) = 0.0;
        
        for n = 1:8
           avail_ptlist(n,1:3) = pt(1:3,n);
           if(pt(4,n) <= 0.0)
               added_points = added_points + 1;
               point_list_for_tet_building(added_points,1:3) = pt(1:3,n);
               point_to_stored_vertex_map(added_points) = n;
           end
        end
        for n = 1:ninterp
            added_points = added_points + 1;
            midpoint = 0.5*(pt(1:3,table_v(iplic,1,n))+pt(1:3,table_v(iplic,2,n)));
            point_list_for_tet_building(added_points,1:3) = midpoint;
            avail_ptlist(8+n,1:3) = midpoint;
            point_to_stored_vertex_map(added_points) = 8+n;
        end
        clear TetsUnder;
        if(added_points > 0)
            TetsUnder = delaunay(point_list_for_tet_building(1:added_points,1), ...
            point_list_for_tet_building(1:added_points,2),...
            point_list_for_tet_building(1:added_points,3));
            number_of_tets_under_plane(iplic) = length(TetsUnder(:,1));
            for n = 1:number_of_tets_under_plane(iplic)
                tet_vertices(iplic,n,1:4) = [point_to_stored_vertex_map(TetsUnder(n,1)), ...
                    point_to_stored_vertex_map(TetsUnder(n,2)),...
                    point_to_stored_vertex_map(TetsUnder(n,3)),...
                    point_to_stored_vertex_map(TetsUnder(n,4))];                
            end
            if(iplic == 16)
                for n = 1:number_of_tets_under_plane(iplic)
                    fprintf('%f %f %f %f \n',...
                    point_to_stored_vertex_map(TetsUnder(n,1)), ...
                    point_to_stored_vertex_map(TetsUnder(n,2)),...
                    point_to_stored_vertex_map(TetsUnder(n,3)),...
                    point_to_stored_vertex_map(TetsUnder(n,4)));
                    figure(n)
                    tetramesh(TetsUnder(n,:),point_list_for_tet_building);
                    hold on;
                    for p = 1:4
                        ind = point_to_stored_vertex_map(TetsUnder(n,p));
                        scatter3(avail_ptlist(ind,1),...
                            avail_ptlist(ind,2),...
                            avail_ptlist(ind,3),40,'MarkerFaceColor','k');
                    end
                    axis([0 1 0 1 0 1]);
                    xlabel('x');
                    ylabel('y');
                    zlabel('z');
                    camup([0 1 0])
                end
                for n = 1:8+ninterp
                   fprintf('Pt #%d:  %f %f %f \n',n,avail_ptlist(n,1),avail_ptlist(n,2),avail_ptlist(n,3)) ;
                end

            end
        end
        
        
        % Collect interface polygon points and positive
        % cube vertex points
        added_points = 0;
        clear point_list_for_tet_building;
        clear point_to_stored_vertex_map;
        point_list_for_tet_building(:,1:3) = 0.0;
        point_to_stored_vertex_map(:) = 0;
        for n = 1:8
           if(pt(4,n) > 0.0)
               added_points = added_points + 1;
               point_list_for_tet_building(added_points,1:3) = pt(1:3,n);
               point_to_stored_vertex_map(n) = n;
           end
        end
        for n = 1:ninterp
            added_points = added_points + 1;
            midpoint = 0.5*(pt(1:3,table_v(iplic,1,n))+pt(1:3,table_v(iplic,2,n)));
            point_list_for_tet_building(added_points,1:3) = midpoint;
            point_to_stored_vertex_map(added_points) = 8+n;
        end
        clear TetsOver;
        if(added_points > 0)
            TetsOver = delaunay(point_list_for_tet_building(1:added_points,1), ...
            point_list_for_tet_building(1:added_points,2),...
            point_list_for_tet_building(1:added_points,3));
            total_number_of_tets(iplic) = max(number_of_tets_under_plane(iplic),0) ...
                +length(TetsOver(:,1));
            for n = 1:length(TetsOver(:,1))
                tet_vertices(iplic,n+max(number_of_tets_under_plane(iplic),0),1:4) = ...
                   [point_to_stored_vertex_map(TetsOver(n,1)), ...
                    point_to_stored_vertex_map(TetsOver(n,2)),...
                    point_to_stored_vertex_map(TetsOver(n,3)),...
                    point_to_stored_vertex_map(TetsOver(n,4))];
            end
        else
            total_number_of_tets(iplic) = number_of_tets_under_plane(iplic);
        end
    end

    % Turn off visualization
    visualize=0;
        
end
 % Hard set this because it is a valid value, just has 0 tets
 % under plane. Therefore should be 0 and not -1.
 number_of_tets_under_plane(256) = 0;
 
 % Now lets write all this out as a C lookup table
 i = 0;
 fprintf('\n constexpr int number_of_tets_under_plane[256] = { \n')
 for ca2 = 1:256
     i = i + 1;
             fprintf('%d, ',number_of_tets_under_plane(ca2));
     if(i >= 16)
        fprintf('\n')
        i = 0;
     end
 end
 fprintf('};');
 
  i = 0;
 fprintf('\n constexpr int total_number_of_tets[256] = { \n')
 for ca2 = 1:256
     i = i + 1;
             fprintf('%d, ',total_number_of_tets(ca2));
     if(i >= 16)
        fprintf('\n')
        i = 0;
     end
 end
 fprintf('};')
 
 fprintf('\n constexpr int tet_vertices[256][17][4] = { \n')
 for ca2 = 1:256
     fprintf('{');
     for s = 1:17
         fprintf('{');
        for v = 1:4
            if(v < 4)
                fprintf('%d, ',tet_vertices(ca2,s,v)-1);
            else
                fprintf('%d ',tet_vertices(ca2,s,v)-1);
            end
        end
        if(s < 17)
            if(mod(s,8) == 0)
                fprintf('}, \n');
            else
                fprintf('}, ');
            end
        else
            fprintf('}');
        end
     end
     fprintf('}, \n');
 end
 fprintf('};\n')
 
 
 

