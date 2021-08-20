//===================================
// Exercise 5: Internal waves
//===================================
//Author: Jochen Kaempf, 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [700,400]; scf(0);

// read input data
rho1=read("rho.dat",-1,101); u1=read("u.dat",-1,101);  
w1=read("w.dat",-1,101); q1=read("q.dat",-1,101);  
[ntot nx] = size(q1); x = (0:5:500)'; z = (0:2:100)'; 
ntot = int(ntot/51); // number of output blocks

for n = 1:60  // animation loop 

time = 10*n; // time in seconds

// grab data blocks
itop = (n-1)*51+1; ibot = itop+50; 
rho = rho1(itop:ibot,1:101)'; u = u1(itop:ibot,1:101)'; 
w = w1(itop:ibot,1:101)'; q = q1(itop:ibot,1:101)'; 

drawlater; clf(); 

// 2d color plot of the density field
Sgrayplot(x,-z,rho,zminmax=[0.0,max(rho)]);
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;500,0];
// contour plot of density field
xset("fpf"," "); col(1:10) = 80; xset("thickness", 2);
contour2d(x,-z,rho,10,col);

xfpoly([148 148 202],[-100 -56 -100]);
xfpoly([0 0 148],[-100 -56 -56]);
xfpoly([0 148 148],[-100 -56 -100]);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;500,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [4,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 100 200 300 400 500], ["0" "100" "200" "300" "400" "500"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]); 

title("Time = "+string(0.1*int(10*time/60))+" min","fontsize",4,'position',[200 -0.5]); // add title

xstring(234, -10,"x (m)");  // add x label
// change font size & set textcolor to white
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(10, -74,"z (m)");  // add z label
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
// end
//end

end // end reference for animation loop


