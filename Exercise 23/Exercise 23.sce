//=======================================
// Exercise 23: Coastal Upwelling in 3D
//=======================================
// Animation of surface distributions of Eulerian concentration & lateral flow field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,500]; scf(0);

// manipulate color map to make the extreme red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

// read input data
c1=read("cS.dat",-1,100); u1=read("uS.dat",-1,100); v1=read("vS.dat",-1,100);
h1=read("h.dat",-1,100);
[ntot nx] = size(c1); x = (1:2:199)'; y = (1:2:50*2-1)';
ntot = int(ntot/50)

for n = 1:ntot // animation loop

nn = n-1; time = nn*6; // time in hrs 

// grab data blosks
itop = (n-1)*50+1; ibot = itop+49; 
cS = c1(itop:ibot,1:100)'; uS = u1(itop:ibot,1:100); vS = v1(itop:ibot,1:100); //grap data block

// interpolation of velocity components onto scalar grid points
um = uS; vm = vS;
for j = 1:50; for k = 2:100; um(j,k) = 0.5*(uS(j,k)+uS(j,k-1)); end; end;
for j = 1:50; um(j,1) = um(j,2); end;
for j = 2:50; for k = 1:100; vm(j,k) = 0.5*(vS(j,k)+vS(j-1,k)); end; end;
for k = 1:100; vm(1,k) = vm(2,k);end;

ua = um(1:4:50,1:4:100); va = vm(1:4:50,1:4:100);

drawlater; clf();

// 2d color plot of Eulerian concentration at the surface
Sgrayplot(x,y,cS,zminmax=[0,0.6]);

// overlay contours
xset("fpf"," "); col(1:11) = 80;  
contour2d(x,y,cS,[0:0.1:1],col);

// overlay flow arrows
champ(x(1:4:100),y(1:4:50),ua',va',1.5);

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

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",4,'position',[100 100]); // add title

xset("color",-1)

for j = 1:50
for k = 1:100
  if h1(j,k) < 10; xfrect(2*(k-1),2*(j-1),2,2); end;
end;
end;

xfrect(0,100,200,2);

xset("color",-1)
//xpoly([0 0],[0 200]);

drawnow;

//if nn < 10 then
//  xs2gif(0,'ex100'+string(nn)+'.gif')
//else
// if nn < 100 then
//    xs2gif(0,'ex10'+string(nn)+'.gif')
//else
//  xs2gif(0,'ex1'+string(nn)+'.gif')
// end
//end



end // end reference for animation loop
