//=============================================
// Exercise 22: Exchange flow through a strait
//=============================================
// Animation of surface and bottom distributions of Eulerian concentrations & flow fields
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [600,600]; scf(0);

// manipulate color map to make the extreme red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

// read input data
C1=read("cS.dat",-1,100); u1=read("uS.dat",-1,100); v1=read("vS.dat",-1,100);
C2=read("cB.dat",-1,100); u2=read("uB.dat",-1,100); v2=read("vB.dat",-1,100);
[ntot nx] = size(C1); x = (1:2:199)'; y = (1:2:99)';
ntot = int(ntot/50);

for n = 1:ntot // animation loop

time = (n-1)*6; // time in hours
nn = n-1;

// grab data blocks
itop = (n-1)*50+1; ibot = itop+49; 
cS = C1(itop:ibot,1:100)'; uS = u1(itop:ibot,1:100); vS = v1(itop:ibot,1:100);
cB = C2(itop:ibot,1:100)'; uB = u2(itop:ibot,1:100); vB = v2(itop:ibot,1:100); 

// subplot 1: surface distribution

// interpolate velocity components onto scalar grid point
um = uS; vm = vS;
for j = 1:50; for k = 2:100; um(j,k) = 0.5*(uS(j,k)+uS(j,k-1)); end; end;
for j = 1:50; um(j,1) = um(j,2); end;
for j = 2:50; for k = 1:100; vm(j,k) = 0.5*(vS(j,k)+vS(j-1,k)); end; end;
for k = 1:100; vm(1,k) = vm(2,k);end;

// elimination of grid points of too low speeds
for j = 1:50; for k = 1:100;
speed(j,k) = sqrt(um(j,k)*um(j,k)+vm(j,k)*vm(j,k));
if speed(j,k) < 0.01; um(j,k) = 0.0; vm(j,k) = 0.0; end;
end; end;

ua = um(1:2:50,1:2:100); va = vm(1:2:50,1:2:100);

drawlater; clf();

subplot(211);

// 2d color plot of surface concentration field
Sgrayplot(x,y,1-cS,,zminmax=[0,1]);

// overlay contour plot 
xset("fpf"," "); col(1:11) = 80;
contour2d(x,y,cS,11,col);

//overlay flow arrows
champ(x(1:2:100),y(1:2:50),ua',va',1);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;200,100];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 20 40 60 80 100], ["0" "20" "40" "60" "80" "100"]);

xset("color",-1)
xfrect(80,30,40,29);
xfrect(80,100,40,30);
xset("color",-1);
xstring(83, 90,"SURFACE");
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[100 100]); // add title
xstring(90,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;
xstring(2,47,"y (cm)");  // add z label
txt=gce(); txt.font_size = 3; txt.font_foreground = -2;

// subplot 2: bottom distribution

// interpolate velocity components onto scalar grid point
um = uB; vm = vB;
for j = 1:50; for k = 2:100; um(j,k) = 0.5*(uB(j,k)+uB(j,k-1)); end; end;
for j = 1:50; um(j,1) = um(j,2); end;
for j = 2:50; for k = 1:100; vm(j,k) = 0.5*(vB(j,k)+vB(j-1,k)); end; end;
for k = 1:100; vm(1,k) = vm(2,k);end;

// elimination of grid points of too low speeds
for j = 1:50; for k = 1:100;
speed(j,k) = sqrt(um(j,k)*um(j,k)+vm(j,k)*vm(j,k));
if speed(j,k) < 0.01; um(j,k) = 0.0; vm(j,k) = 0.0; end;
end; end;

ua = um(1:2:50,1:2:100); va = um(1:2:50,1:2:100);

subplot(212);

// 2d color plot of bottom concentration field
Sgrayplot(x,y,1-cB,,zminmax=[0,1]);

// overlay contour plot 
xset("fpf"," "); col(1:11) = 80;
contour2d(x,y,cB,11,col);

//overlay flow arrows
champ(x(1:2:100),y(1:2:50),ua',va',1);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;200,100];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 40 80 120 160 200], ["0" "40" "80" "120" "160" "200"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 20 40 60 80 100], ["0" "20" "40" "60" "80" "100"]);

xset("color",-1)
xfrect(80,30,40,29);
xfrect(80,100,40,30);
xset("color",-1);
xstring(83, 90,"BOTTOM");
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",3,'position',[100 100]); // add title
xstring(90,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1;
xstring(2,47,"y (cm)");  // add z label
txt=gce(); txt.font_size = 3; txt.font_foreground = -2;

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
