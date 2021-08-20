//==================================
// Bathymetry map for Exercise 22
//==================================

f = gcf(); f.pixmap='on'; 

// reverse grayscale map
map = 1-graycolormap(64); 
f.color_map = map;

f.figure_size = [1000,500]; // set size of graphic window

// read bathymetry data
h1=read("topo.dat",-1,102); 
x = (0:2:202)'; y = (0:2:102)';  
hzero = max(10*int(h1/10),0.0); // show the distribution with 10-m depth increments

// draw 2d grayscale map
Sgrayplot(x,y,hzero',zminmax=[0,140]);

// overlay contours
xset("fpf"," "); b(1:11) = 80;
contour2d(x,y,hzero',[0:10:100],b,strf='211');

// specify graph & axis properties
a = get("current_axes"); 
a.font_size = 3; a.data_bounds = [0,0;200,100];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 20 40 60 80 100], ["0" "20" "40" "60" "80" "100"]);

xstring(90,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(2,47,"y (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

show_pixmap(); 

