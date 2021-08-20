//==================================
// Bathymetry map for Exercise 23
//==================================

f = gcf(); f.figure_size = [1000,500]; scf(0);
// reverse grayscale map
map = 1-graycolormap(64); 
f.color_map = map;

drawlater; clf();
// read bathymetry data
h1=read("topo.dat",-1,102);
x = (0:2:202)'; y = (0:2:102)'; 
hzero = max(10*int(h1/10),0.0); // convert data on the basis of 10-m depth increments

// draw 2d grayscale map
Sgrayplot(x,y,hzero',zminmax=[0,140]);

// overlay contours
b(1:10) = 80; xset("fpf"," ");
contour2d(x,y,hzero',[20:10:110]);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;200,100];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 20 40 60 80 100], ["0" "20" "40" "60" "80" "100"]);

xstring(100,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;
xstring(2,27,"y (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -2;


xset("color",32)
// add land mask
for j = 1:50
for k = 1:100
 if h1(j,k) < 10; xfrect(2*(k-1),2*(j-1),2,2); end;
end;
end;
xfrect(0,100,200,2);

xset("color",-1)

drawnow; 
