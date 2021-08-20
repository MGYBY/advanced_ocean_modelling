//===========================================
// Exercise 21: Eddy formation in a Strait
//===========================================
// Animation of Eulerian concentration & surface flow fields
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [1000,500]; scf(0);

// manipulate color map to make the dark red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

// read input data
c1=read("cS.dat",-1,100); u1=read("uS.dat",-1,100); v1=read("vS.dat",-1,100);
[ntot nx] = size(c1); x = (1.8:3.6:100*3.6-1.8)'; y = (1.8:3.6:50*3.6-1.8)'; 
ntot = int(ntot/50)

for n = 1:ntot // animation loop

time = (n-1)*3; // time in hours
nn = n-1;

// grab data blocks
itop = (n-1)*50+1; ibot = itop+49; 
cS = c1(itop:ibot,1:100)'; uS = u1(itop:ibot,1:100); vS = v1(itop:ibot,1:100);

// interpolation of velocity components onto scalar grid points
um = uS; vm = vS;
for j = 1:50; for k = 2:100; um(j,k) = 0.5*(uS(j,k)+uS(j,k-1)); end; end;
for j = 1:50; um(j,1) = um(j,2); end;
for j = 2:50; for k = 1:100; vm(j,k) = 0.5*(vS(j,k)+vS(j-1,k)); end; end;
for k = 1:100; vm(1,k) = vm(2,k);end;

// elimination of grid points of small speeds
for j = 1:50; for k = 1:100;
speed(j,k) = sqrt(um(j,k)*um(j,k)+vm(j,k)*vm(j,k));
if speed(j,k) < 0.01; um(j,k) = 0.0; vm(j,k) = 0.0; end;
end; end;

drawlater; clf;

// 2d colour plot of concentration field
Sgrayplot(x,y,cS,zminmax=[0,0.5]);

// overlay concentration contours
xset("fpf"," "); col(1:11) = 80;  
contour2d(x,y,cS,11,col);

ua = um(1:2:50,1:2:100); va = um(1:2:50,1:2:100);

// overlay flow arrows
champ(x(1:2:100),y(1:2:50),ua',va',1);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;360,180];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [1,1];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 60 120 180 240 300 360], ["0" "60" "120" "180" "240" "300" "360"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 60 120 180], ["0" "60" "120" "180"]);

xset("color",-1)
xfrect(72,28,288,27);
xfrect(0,180,33,108);

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",4,'position',[180 180]); // add title
xstring(135,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(1,100,"y (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;

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
