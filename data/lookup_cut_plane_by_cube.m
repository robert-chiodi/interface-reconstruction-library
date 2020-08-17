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
    pt(1,2) = x(2); pt(2,2) = y(2); pt(3,2) = z(1); 
    pt(1,3) = x(2); pt(2,3) = y(2); pt(3,3) = z(2);
    pt(1,4) = x(2); pt(2,4) = y(1); pt(3,4) = z(2);
    pt(1,5) = x(1); pt(2,5) = y(1); pt(3,5) = z(1);
    pt(1,6) = x(1); pt(2,6) = y(2); pt(3,6) = z(1);
    pt(1,7) = x(1); pt(2,7) = y(2); pt(3,7) = z(2);
    pt(1,8) = x(1); pt(2,8) = y(1); pt(3,8) = z(2);

% Data to tabulate
table_done=zeros(1,256);
table_ninterp(1:256)=-1;
table_v(1:256,1:2,1:6)=0;

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
        
        % Test bogus case
        if iplic==222
            visualize=1;
        end
        
        % Display the problem
        for n=1:8
            vertices(n,1:3)=pt(1:3,n);
        end
        faces=[1 2 3 4; 5 6 7 8; 1 5 6 2; 4 3 7 8; 1 4 8 5; 2 6 7 3];
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
        
    end

    % Turn off visualization
    visualize=0;
        
end

    % Now lets write all this out as a C lookup table
 i = 0;
 fprintf('\n constexpr int number_of_plane_intersections_with_cuboid[256] = { \n')
 for ca2 = 1:256
     i = i + 1;
             fprintf('%d, ',table_ninterp(ca2));
     if(i >= 16)
        fprintf('\n')
        i = 0;
     end
 end
 fprintf('};')
 
 fprintf('\n constexpr int cuboid_cut_vertices[256][2][6] = { \n')
 for ca2 = 1:256
     fprintf('{');
     for s = 1:2
         fprintf('{');
        for v = 1:6
            if(v < 6)
                fprintf('%d, ',table_v(ca2,s,v)-1);
            else
                fprintf('%d ',table_v(ca2,s,v)-1);
            end
        end
        if(s == 1)
            fprintf('}, ');
        end
        if(s == 2)
            fprintf('}');
        end
     end
     fprintf('}, \n');
 end
 fprintf('};\n')
