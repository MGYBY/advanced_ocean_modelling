//============================
// Exercise 19: Ekman pumping
//============================
// Animation of sea level & density field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,600]; scf(0);

// read input data
ep=read("eta.dat",-1,101);
rho1=read("rho.dat",-1,101); 

[ntot nx] = size(ep); x = (0:5:500)'; z = (10:20:510)';

for n = 1:ntot // animation loop 

//*********************************
// Top graph: sea-level elevation
//*********************************

time = (n-1)*12; 

// grab data blocks
itop = (n-1)*26+1; ibot = itop+25; 
rho = rho1(itop:ibot,1:101)';

drawlater; clf;

subplot(211)
plot2d(x,100*ep(n,:),3);
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;500,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 100 200 300 400 500], ["0" "100" "200" "300" "400" "500"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -50 0], ["-100" "-50" "0"]); 
a.axes_bounds = [0. 0. 1. 0.2];
xstring(345, -94,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(3, -55,"eta (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

subplot(212)
// 2d color plot of density field
Sgrayplot(x,-z,7-rho,zminmax=[0,8]);
// overlay density contours
xset("fpf"," "); col(1:15) = 80; 
contour2d(x,-z,rho,[0.:0.5:7],col);
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-500;500,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 100 200 300 400 500], ["0" "100" "200" "300" "400" "500"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-500 -400 -300 -200 -100 0], ["-500" "-400" "-300" "-200" "-100" "0"]); 
a.axes_bounds = [0. 0.15 1. 0.85];

xstring(345, -490,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(3, -260,"z (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[410 -10]);

drawnow;

// save frames as sequential GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
// if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//else
//  xs2gif(0,'ex1'+string(n)+'.gif')
// end
//end

end // end reference for animation loop
