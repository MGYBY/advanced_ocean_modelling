//============================================
// Exercise 20: Geostrophic adjustment in 3d
//============================================
// Animation of density & transverse flow field in cross section
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,300]; scf(0);

// read input data
rho1=read("rhoV.dat",-1,25); v1=read("vV.dat",-1,25);
[ntot nx] = size(rho1); x = (1:2:49)'; z = (10:20:490)'; 
ntot = int(ntot/25);

for n = 1:ntot // animation loop

nn = n-1; time = (n-1); // time in hours

//grab data blocks 
itop = (n-1)*25+1; ibot = itop+24; 
rho = rho1(itop:ibot,1:25)'; v = v1(itop:ibot,1:25)'; 

drawlater; clf();

// draw 2d color plot of density field
Sgrayplot(x,-z,rho,zminmax=[-0.1,0.01]);

// draw contours of v distinguishing positive and negative values
xset("fpf"," "); col(1:10) = 80;  
contour2d(x,-z,v,[0.02:0.04:0.4],col);

xset("fpf"," "); col(1:10) = 80;  xset("line style",2)
contour2d(x,-z,-v,[0.02:0.04:0.4],col);
xset("line style",1)

xstring(24, -480,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(3, -260,"z (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-500;50,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-500 -400 -300 -200 -100 0], ["-500" "-400" "-300" "-200" "-100" "0"]); 

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",4,'position',[25 0]); // add title

drawnow;

// save frames as sequential GIF files (optional)
//if nn < 10 then
//  xs2gif(0,'ex100'+string(nn)+'.gif')
//else
// if nn < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//else
//  xs2gif(0,'ex1'+string(n)+'.gif')
// end
//end

end // end reference for animation loop
