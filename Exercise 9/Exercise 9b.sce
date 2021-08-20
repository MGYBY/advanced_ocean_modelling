//======================================
// Exercise 9: Convective Entrainment
//======================================
// This script produces an animation of float movements
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.figure_size = [1000,400]; scf(0);

// read input data
trx1=read("TRx.dat",-1,3000); // x-coordinates of floats
try1=read("TRz.dat",-1,3000); // z-coordinates of floats
rho1=read("rho.dat",-1,201); // density field
[ntot nx] = size(rho1); x = (0:5:1000)'; z = (0:5:100)';
[ntot ntra] =size(trx1);
time = 0.0;

for n = 1:ntot // animation loop

time = n*3; // time in minutes

// grab data block 
itop = (n-1)*21+1; ibot = itop+20; 
rho = rho1(itop:ibot,1:201)'; 
 // grab data line
trax = trx1(n,1:ntra); tray = try1(n,1:ntra); 

drawlater; clf(); 

// draw floats
plot2d(trax,-tray+2,-3);//
m = gce(); m.children.mark_size = 0.1;

xset("fpf"," "); col(1:10) = 4;
contour2d(x,-z,rho,10,col);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;1000,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

title("Time = "+string(0.1*int(10*time))+" min","fontsize",4,'position',[400 0]); // add title

xstring(470, -15,"x (m)");  // add x label
// change font size & set textcolor to black
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(5, -54,"z (m)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

drawnow;

// save frames as sequential GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end // end reference for animation loop
