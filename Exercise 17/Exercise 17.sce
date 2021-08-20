//==============================================
// Exercise 17: Tidal Mixing Fronts
//==============================================
// Animation of density field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1200,400]; scf(0);

// read input data
rho1=read("rho.dat",-1,101); 
[ntot nx] = size(rho1); 
x = (0:0.1:10)'; z = (1:2:101)';
ntot = int(ntot/51);

for n = 1:ntot // animation loop

time = (n-1)*3; // time in hours

// grab data blocks
itop = (n-1)*51+1; ibot = itop+50; 
rho = rho1(itop:ibot,1:101)'; 

// remove artificak gradients near seafloor
for k = 1:101
for i = 11:21
if rho(k,i) == 0.0; rho(k,i) = rho(k,i-1); end;
end;
end;

drawlater; clf();

// color plot of density field
Sgrayplot(x,-z,1.0-rho,zminmax=[0,1]); //swapping colors
colorbar(0,1);
cc = gce(); cc.parent.font_size = 2;

// overlay contour plot of density field
xset("fpf"," "); col(1:16) = 80; 
contour2d(x,-z,rho,[0.:0.1:1],col);
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;10,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 2 4 6 8 10], ["0" "2" "4" "6" "8" "10"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

// bathymetry mask
xfpoly([0,0,10,10],[-100,-100,-48,-100]);

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[5 0]); // add title

xstring(6.5, -97,"x (m)");  // add x label
// change font size & set textcolor to black
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(0.3, -55,"z (m)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;

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
