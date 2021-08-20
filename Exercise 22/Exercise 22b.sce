//=============================================
// Exercise 22: Exchange flow through a strait
//=============================================
// Animation of density distribution and tranverse flow in a vertical transect cutting across the middle of the strait
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [600,600]; scf(0);

// manipulate color map to make the extreme red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

f.figure_size = [800,300]; // set size of graphic window

// read input data
rho1=read("rhoV.dat",-1,50); u1=read("uV.dat",-1,50); h1=read("topo.dat",-1,102);
depth(1:50) = 10*floor(h1(2:51,51)/10); // create coarse bathymetry mask
[ntot ny] = size(rho1); y = (1:2:99)'; z = (5:10:105)'; 
ntot = int(ntot/11)

for n = 1:ntot // animation loop

nn = n-1; time = (n-1)*6; // time in hours

// grab data blocks
itop = (n-1)*11+1; ibot = itop+10; 
rho = rho1(itop:ibot,1:50)'; u = u1(itop:ibot,1:50)';

drawlater; clf();

// draw 2d color plot of density field
Sgrayplot(y,-z,rho,zminmax=[-0.05,0.55]);

// overlay contours of v
xset("fpf"," "); col(1:10) = 80;
contour2d(y,-z,-u,[0.05:0.05:0.5],col); // draw negative values with dashed lines
xset("fpf"," "); col(1:10) = 80; xset("line style",1) 
contour2d(y,-z,u,[0.05:0.05:0.5],col); // draw positive values with full lines

// draw bathymetry
xset("color",33)
for k = 1:50
  xx = real(k-1)*2;
  yy = -depth(k);
  width = 2;
  height = 110.-yy;
  xfrect(xx,yy,width,height);
end;
xset("color",0)

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-100;100,0];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 20 40 60 80 100], ["0" "20" "40" "60" "80" "100"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]);

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[50 0]); // add title
xstring(48,-96,"y (m)");  // add x label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;
xstring(1,-54,"z (cm)");  // add z label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;

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


