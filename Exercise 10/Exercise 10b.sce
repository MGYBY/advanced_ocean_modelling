//==============================================
// Exercise 10: Slope convection near the shore
//==============================================
// Animation of the evolution of the Eulerian concentration field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0);  

// read input data
c1=read("c.dat",-1,201); 
[ntot nx] = size(c1); x = (0:5:1000)'; z = (2.5:5:102.5)';
ntot = int(ntot/21);

for n = 1:120 // animation loop 

time = 6*n; // time in minutes
// grab data blocks
itop = (n-1)*21+1; ibot = itop+20; 
c = c1(itop:ibot,1:201)'; 

drawlater; clf(); 

// 2d colour plot of density field
Sgrayplot(x,-z,c,zminmax=[0.0,1]);
// contour plot of density field
xset("fpf"," "); col(1:20) = 80;
// selected contours
d = [0.05:0.05:1];
contour2d(x,-z,c,d,col);
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;1000,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

// draw bathymetry mask
xfpoly([0 0 260 260],[-18 -100 -100 -18]);
xfpoly([260 260 510],[-18 -100 -100]);
xfpoly([0 0 10 10],[0 -18 -18 0]);
xfpoly([990 990 1000 1000],[0 -100 -100 0]);
xfpoly([0 0 1000 1000],[-98 -100 -100 -98]);

title("Time = "+string(0.1*int(10*time))+" min","fontsize",4,'position',[400 0]); // add title

xstring(390, -97,"x (m)");  // add x label
// change font size & set textcolor to black
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
