//==============================
// Exercise 8: Free convection
//==============================
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0);

// read input data (only density here)
rho1=read("rho.dat",-1,201); 
[ntot nx] = size(rho1); x = (0:5:1000)'; z = (0:5:100)'; 
ntot = int(ntot/21);

for n = 1:ntot

time = 5*n; // time in minutes

// grab data block (only density here)
itop = (n-1)*21+1; ibot = itop+20; 
rho = rho1(itop:ibot,1:201)'; 

drawlater; clf(); 
// 2d color plot of the density field
Sgrayplot(x,-z,rho,zminmax=[0,0.02]);

// 2d contour plot of the density field
xset("fpf"," "); col(1:5) = 80;
contour2d(x,-z,rho,5,col);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;1000,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

title("Time = "+string(0.1*int(10*time))+" min","fontsize",4,'position',[400 0]); // add title

xstring(470, -97,"x (m)");  // add x label
// change font size & set textcolor to black
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(10, -54,"z (m)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

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
