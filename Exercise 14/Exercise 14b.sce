//=================================
// Exercise 14: Positive estuaries
//=================================
// Animation of age distribution and flow field
//Author: Jochen Kaempf, 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1200,400]; scf(0);

// read input data
f1=read("age.dat",-1,101); u1=read("uM.dat",-1,101); w1=read("wM.dat",-1,101);
// read bathymetry data
dep=read("h.dat",1,101);

[ntot nx] = size(f1); x = (0:2:200)'; z = (0.5:1:20.5)'; 
ntot = int(ntot/21);

for n = 1:ntot // snapshot output
  
time = n/2; // time in hours

// grab data blocks
itop = (n-1)*21+1; ibot = itop+20; 
f = f1(itop:ibot,1:101); u = u1(itop:ibot,1:101); w = w1(itop:ibot,1:101); 

// manipulate age in dry grid cells for visualisation of bottom topography
for k = 1:101
 depth(k) = 5.0 + 15.0*(k-1)/100.0;
 nb(k) = floor(depth(k)/1.0);
 f(nb(k)+1,k)= f(nb(k),k);
 nb(k) = min(nb(k),19);
 f(nb(k)+2,k)= f(nb(k),k);
end;

drawlater; clf;

// 2d color plot of agea distribution
Sgrayplot(x,-z,f',zminmax=[0 4]);
colorbar(0,4);
cc = gce(); cc.parent.font_size = 3;
cc.parent.title.text = "Age (days)",cc.parent.title.font_size = 3;
// contour plot of age distribution
xset("fpf"," "); col(1:8) = 80; 
contour2d(x,-z,f',[0.5:0.5:4],col);

ua = u(1:2:21,1:5:101); wa = w(1:2:21,1:5:101);
// draw velocity vector field
champ(x(1:5:101),-z(1:2:21),ua',wa',2.0); 
// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-20;200,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-20 -15 -10 -5 0], ["-20" "-15" "-10" "-5" "0"]); 

title("Time = "+string(0.01*int(100*time))+" days","fontsize",3,'position',[95 0]); // add title

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

