//=================================
// Exercise 14: Positive estuaries
//=================================
// Animation of density field and currents
//Author: Jochen Kaempf, 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,400]; scf(0);

// read input data
rho1=read("rho.dat",-1,101); u1=read("u.dat",-1,101); w1=read("w.dat",-1,101);
// read bathymetry data
dep=read("h.dat",1,101);

[ntot nx] = size(rho1); x = (0:2:200)'; z = (0.5:1:20.5)';
ntot = int(ntot/21);

for n = 1:ntot // animation loop

time = (n-1)*2; // time in hours

// grab data blocks
itop = (n-1)*21+1; ibot = itop+20; 
rho = rho1(itop:ibot,1:101); u = u1(itop:ibot,1:101); w = w1(itop:ibot,1:101); 

// manipulate density in dry grid cells for visualisation of bottom topography
for k = 1:101
 depth(k) = 5.0 + 15.0*(k-1)/100.0;
 nb(k) = floor(depth(k)/1.0);
 rho(nb(k)+1,k)= rho(nb(k),k);
 nb(k) = min(nb(k),19);
 rho(nb(k)+2,k)= rho(nb(k),k);
end;

drawlater; clf;

// 2d color graph of density field
Sgrayplot(x,-z,rho',zminmax=[0.0,max(rho)]);

// draw density contours
xset("fpf"," "); col(1:15) = 80;
contour2d(x,-z,rho',[0:2:28],col);

ua = u(1:2:21,1:5:101); wa = w(1:2:21,1:5:101);

// draw velocity vector field
champ(x(1:5:101),-z(1:2:21),ua',wa',2.0); 
// specify graph & axis properties
a = gca(); 
a.font_size = 3; a.data_bounds = [0,-20;200,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-20 -15 -10 -5 0], ["-20" "-15" "-10" "-5" "0"]); 

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[95 0]); // add title

xfpoly([0 0 200],[-5 -20 -20]);

xstring(95, -19.5,"x (m)");  // add x label
// change font size & set textcolor to black
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(2, -11,"z (m)");  // add z label
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
