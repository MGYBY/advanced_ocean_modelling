//=====================================
// Display bathymetry for Exercise 21
//=====================================

f = gcf(); f.color_map = 1-graycolormap(64);
f.figure_size = [800,600]; // set size of graphic window

// read input data
h1=read("topo.dat",-1,102); // read input data
x = (0:1*3.6:101*3.6)'; y = (0:1*3.6:51*3.6)'; // location vectors  
hzero = max(h1,0.0);

drawlater; clf();

Sgrayplot(x,y,hzero',zminmax=[0,300]);
col(1:10) = 200; xset("fpf"," ");
contour2d(x,y,hzero',[20:20:200],col);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;360,180];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 60 120 180 240 300 360], ["0" "60" "120" "180" "240" "300" "360"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 60 120 180], ["0" "60" "120" "180"]);

xstring(135,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(2,83,"y (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
drawnow(); 


