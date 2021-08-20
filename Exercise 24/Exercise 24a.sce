//=======================================
// Exercise 24: The abyssal circulation
//=======================================
// Animation of bottom distributions of density & lateral flow field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [700,700]; scf(0);

// manipulate color map to make the dark red a bit lighter
map = jetcolormap(64); 
ic = 57; for i = ic+1:64; map(i,1:3) = map(ic,1:3); end;
f.color_map = map;

// read input data
rho2=read("rhoB.dat",-1,51); u2=read("uB.dat",-1,51); v2=read("vB.dat",-1,51);
[ntot nx] = size(rho2); x = (10:20:1010)'; y = (10:20:1010)';
ntot = int(ntot/51);

for n = 1:ntot // animation loop

nn = n-1; time = nn;

// grab data blocks
itop = (n-1)*51+1; ibot = itop+50; 
rhoB = rho2(itop:ibot,1:51)'; uB = u2(itop:ibot,1:51); vB = v2(itop:ibot,1:51);

// interpolation of velocity components onto scalar grid points
um = uB; vm = vB;
for j = 1:51; for k = 2:51; um(j,k) = 0.5*(uB(j,k)+uB(j,k-1)); end; end;
for j = 1:51; um(j,1) = um(j,2); end;
for j = 2:51; for k = 1:51; vm(j,k) = 0.5*(vB(j,k)+vB(j-1,k)); end; end;
for k = 1:51; vm(1,k) = vm(2,k);end;

ua = um(1:2:51,1:2:51); va = vm(1:2:51,1:2:51);

// elimination of grid cells with too low speeds
for j = 1:25; for k = 1:25;
speed(j,k) = sqrt(ua(j,k)*ua(j,k)+va(j,k)*va(j,k));
if speed(j,k) < 0.01; ua(j,k) = 0.0; va(j,k) = 0.0; end;
end; end;

drawlater; clf();

// 2d color map of density field
Sgrayplot(x,y,rhoB,zminmax=[0,2]);

// overlay contours
xset("fpf"," "); col(1:11) = 80;  
contour2d(x,y,rhoB,11,col);
 
// overlay velocity arrows
champ(x(1:2:51),y(1:2:51),ua',va',1);

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0;1000,1000];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);

xstring(450,-62,"x (m)");  // add x label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1; txt.clip_state = "off";
xstring(-109,470,"y (cm)");  // add z label
txt=gce(); txt.font_size = 3; txt.font_foreground = -1; txt.clip_state = "off";

title("Time = "+string(0.01*int(100*time))+" days","fontsize",3,'position',[600 1000]); // add title

// overlay land mask
xset("color",0);
xfrect(0,1000,400,200);
xfrect(600,1000,400,200);
xfrect(0,1000,1000,10);
xfrect(0,10,1000,10);
xfrect(0,1000,10,1000);
xfrect(990,1000,10,1000);

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
