//=============================================
// Exercise 25: Simulation of an El-Nino Event
//=============================================
// Animation of distributions of Eulerian concentration at 100 m & lateral flow field at the surface
// Author: Jochen Kaempf, March 2015 (update)

f = gcf();f.figure_size = [700,700]; scf(0);

// manipulate color map to make the dark red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

// read input data
ci=read("c3.dat",-1,101); ui=read("uS.dat",-1,101); vi=read("vS.dat",-1,101);
[ntot nx] = size(ci); x = (0:20:2000)'; y = (-1000:20:1000)'; 
ntot = int(ntot/101)

for n = 1:ntot // animation loop

nn = n-1; time = nn/2; // time in days

// grab data blocks
itop = (n-1)*101+1; ibot = itop+100; 
c = ci(itop:ibot,1:101)'; u = ui(itop:ibot,1:101); v = vi(itop:ibot,1:101);

// interpolate velocity components onto scalar grid point
um = u; vm = v;
for j = 1:101; for k = 2:101; um(j,k) = 0.5*(u(j,k)+u(j,k-1)); end; end;
for j = 1:101; um(j,1) = um(j,2); end;
for j = 2:101; for k = 1:101; vm(j,k) = 0.5*(v(j,k)+v(j-1,k)); end; end;
for k = 1:101; vm(1,k) = vm(2,k);end;

ua = um(1:4:101,1:4:101); va = vm(1:4:101,1:4:101);

// eliminate grid cells of too low speed
for j = 1:25; for k = 1:25;
speed(j,k) = sqrt(ua(j,k)*ua(j,k)+va(j,k)*va(j,k));
if speed(j,k) < 0.02; ua(j,k) = 0.0; va(j,k) = 0.0; end;
end; end;

drawlater; clf();

// 2d color plot of Eulerian concentration
Sgrayplot(x,y,c,zminmax=[0.5,1.02*max(c)]);

// overlay contours
xset("fpf"," "); col(1:15) = 80;
contour2d(x,y,c,15,col);

//xset("color",80);
equator(1:nx) = 0.0;
plot2d(x,equator,80) 
// overlay flow vectors
champ(x(1:4:101),y(1:4:101),ua',va',1);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,-1000;2000,1000];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [4,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 500 1000 1500 2000], ["0" "500" "1000" "1500" "2000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [-1000 -500 0 500 1000], ["-1000" "-500" "0" "500" "1000"]);

xset("color",0);
xfrect(0,1000,30,2000);
xfrect(2000-30,1000,30,2000);

xstring(1150,-1110,"x(km)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1; txt.clip_state = "off";
xstring(-200,40,"y(km)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1; txt.clip_state = "off";

title("Time = "+string(0.01*int(100*time))+" days","fontsize",4,'position',[1000 1000]); // add title

drawnow;

// save frames as sequential GIF files (optional)
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
