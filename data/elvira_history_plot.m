%  This file is part of the Interface Reconstruction Library (IRL),
%  a library for interface reconstruction and computational geometry operations.
%
%  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
%
%  This Source Code Form is subject to the terms of the Mozilla Public
%  License, v. 2.0. If a copy of the MPL was not distributed with this
%  file, You can obtain one at https://mozilla.org/MPL/2.0/.

clear all; close all;
x0 = 0.0;
y0 = 0.0;
z0 = 0.0;
meshsize = 1.0;%/32.0;
run('reconstruction_history.m');
number_of_saved_reconstructions = length(bestPlane(1,1,1,1,:));

for iteration = 1:number_of_saved_reconstructions
   figure(iteration)
   hold on;
   for j = -1 : 1
       for i = -1 : 1
           cell = (j+1)*3+(i+1)+1;
           pt(1,1:3) = [x0-0.5*meshsize,y0-0.5*meshsize,0.0] + [i*meshsize,j*meshsize,0.0];
           pt(2,1:3) = [x0+0.5*meshsize,y0-0.5*meshsize,0.0] + [i*meshsize,j*meshsize,0.0];
           pt(3,1:3) = [x0+0.5*meshsize,y0+0.5*meshsize,0.0] + [i*meshsize,j*meshsize,0.0];
           pt(4,1:3) = [x0-0.5*meshsize,y0+0.5*meshsize,0.0] + [i*meshsize,j*meshsize,0.0];
           handle = fill3(pt(1:4,1),pt(1:4,2),pt(1:4,3),'b');
           alpha(handle,liq_VF(cell));
       end
   end
   axis([x0-1.5*meshsize x0+1.5*meshsize y0-1.5*meshsize y0+1.5*meshsize z0-1.5*meshsize z0+1.5*meshsize]);
   %axis([-0.5*meshsize 0.5*meshsize -0.5*meshsize 0.5*meshsize -0.5*meshsize 0.5*meshsize]);
   if(cell == 9) % plot 2D lines
    lw = 3;
    hline(y0+0.5*meshsize,lw);
    hline(y0+1.5*meshsize,lw);
    hline(y0-0.5*meshsize,lw);
    hline(y0-1.5*meshsize,lw);
    vline(x0+0.5*meshsize,lw);
    vline(x0+1.5*meshsize,lw);
    vline(x0-0.5*meshsize,lw);
    vline(x0-1.5*meshsize,lw);
   end
   for cell = 1:length(liq_VF(:))
       for plane = 1:bestPlaneNPlane(iteration, cell)
           if(bestPlaneNvert(cell,plane,iteration) == 0)
              continue; 
           end
            fill3(bestPlane(1,1:bestPlaneNvert(cell,plane,iteration),cell,plane,iteration),...
            bestPlane(2,1:bestPlaneNvert(cell,plane,iteration),cell,plane,iteration),...
            bestPlane(3,1:1:bestPlaneNvert(cell,plane,iteration),cell,plane,iteration),'k','LineWidth',4);
       end
   end
   hold off;
   view(2)
   if(cell == 9)
    view(2);    
    xlabel('x');
    ylabel('y');
   end
   if(cell == 27)
    view(3);    
    xlabel('x');
    ylabel('y');
    zlabel('z');
   end
   xt = get(gca, 'XTick');
   set(gca, 'FontSize', 1);
   yt = get(gca, 'YTick');
   set(gca, 'FontSize', 1)
   prefix = 'plots/R2P_Accepted_Step_';
   iter_string = sprintf('%04d',iteration);
   pname = sprintf(strcat(prefix,iter_string));
   print(pname,'-dpng')
    
end
