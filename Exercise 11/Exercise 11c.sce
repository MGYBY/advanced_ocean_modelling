//============================================
// Exercise 11: Double-diffusive instability
//============================================
// Animation of density field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0);  

// read input data
rho1=read("rho.dat",-1,201); 
[ntot nx] = size(rho1); x = (0:1:200)'; z = (0.5:1:20.5)';
ntot = int(ntot/21);

for n = 1:20 // animation loop

time = n; // time in minutes

// grab data blocks
itop = (n-1)*21+1; ibot = itop+20; 
rho = rho1(itop:ibot,1:201)'; 

drawlater; clf();

// 2d colour plot of density field
Sgrayplot(x,-z,rho,zminmax=[0.0,max(rho)]);

// contour plot of density field
xset("fpf"," "); col(1:10) = 80;
contour2d(x,-z,rho,10,col);
 
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-20;200,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-20 -16 -12 -8 -4 0], ["-20" "-16" "-12" "-8" "-4" "0"]); 

title("Time = "+string(0.1*int(10*time))+" min","fontsize",3,'position',[90 0]); // add title

xstring(96, -22,"x (m)");  // add x label
b = gce(); b.clip_state = "off";
// change font size & set textcolor to black
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;
xstring(-11, -11,"z (m)");  // add z label
b = gce(); b.clip_state = "off";
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;

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
