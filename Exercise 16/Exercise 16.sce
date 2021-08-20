//=====================================
// Exercise 16: Geostrophic adjustment
//=====================================
// Animation of density field and flow structure
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0);

// read input data
rho1=read("rho.dat",-1,101); u1=read("u.dat",-1,101); v1=read("v.dat",-1,101);
[ntot nx] = size(rho1); x = (0:0.5:50)'; z = (0:10:500)'; 
ntot = int(ntot/51);

for n = 1:ntot // animation loop

time = n; // time in hours

// grab data blocks
itop = (n-1)*51+1; ibot = itop+50; 
rho = rho1(itop:ibot,1:101)'; u = u1(itop:ibot,1:101); v = v1(itop:ibot,1:101); 

drawlater; clf;
// 2d color plot of density field 
Sgrayplot(x,-z,rho,zminmax=[-0.11,0]);

// contour plots of v-velocity component (distinguishing between positive and negative values)
xset("fpf"," "); col1(1:10) = 90; xset("line style",1); xset("thickness",2);
contour2d(x,-z,v',[0.05:0.05:0.5],col1);

xset("fpf"," "); col2(1:11) = 90; xset("line style",2); xset("thickness",2);
contour2d(x,-z,-v',[0.0:0.05:0.5],col2);
 // specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-500;50,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [4,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-500 -400 -300 -200 -100 0], ["-500" "-400" "-300" "-200" "-100" "0"]); 

title("Time = "+string(0.1*int(10*time))+" hrs","fontsize",4,'position',[22 0]); // add title

xstring(24, -480,"x (m)");  // add x label
// change font size & set textcolor to black
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(1, -275,"z (m)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;

drawnow;

// save frames as sequential GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
// if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
// else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end // end reference for animation loop
