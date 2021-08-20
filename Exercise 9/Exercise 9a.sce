//====================================
// Exercise 9: Convective entrainment
//====================================
// Animation of Eulerian tracer concentration
// Author: Jochen Kaempf, March 2015 (update)
f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0); 

// read input data 
rho1=read("rho.dat",-1,201); c1=read("c.dat",-1,201); 
[ntot nx] = size(rho1); x = (0:5:1000)'; z = (0:5:100)'; 
ntot = int(ntot/21);

for n = 1:ntot

time = 3*n; // time in minutes
// grab data blocks
itop = (n-1)*21+1; ibot = itop+20; 
rho = rho1(itop:ibot,1:201)'; c = c1(itop:ibot,1:201);

drawlater; clf();

// 2d color plot of the density field
Sgrayplot(x,-z,c',zminmax=[0,1]);

// 2d contour plot of the density field
xset("fpf"," "); col(1:10) = 80;
contour2d(x,-z,rho,10,col);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;1000,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

title("Time = "+string(0.1*int(10*time))+" min","fontsize",4,'position',[400 0]); // add title

xstring(470, -97,"x (m)");  // add x label
// change font size & set textcolor to white
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(10, -54,"z (m)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;

drawnow;

// save frames as sequential GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
// else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end // end reference for animation loop
