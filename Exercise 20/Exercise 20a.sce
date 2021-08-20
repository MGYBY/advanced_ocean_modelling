//============================================
// Exercise 20: Geostrophic adjustment in 3d
//============================================
// Animation of surface density & flow field
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = jetcolormap(64); f.figure_size = [600,600]; scf(0);

// read input data
rho1=read("rhoS.dat",-1,25); u1=read("uS.dat",-1,25); v1=read("vS.dat",-1,25);
[ntot nx] = size(rho1); x = (1:2:49)'; y = (1:2:49)';
ntot = int(ntot/25);

for n = 1:ntot //animation loop

time = (n-1); // time in hours
nn = n-1;

// grab data blocks
itop = (n-1)*25+1; ibot = itop+24; 
rho = rho1(itop:ibot,1:25)'; u = u1(itop:ibot,1:25)'; v = v1(itop:ibot,1:25)';

// interpolate onto scalar grid point
um = u; vm = v;
for j = 1:25; for k = 2:25; um(j,k) = 0.5*(u(j,k)+u(j,k-1)); end; end;
for j = 1:25; um(j,1) = um(j,2); end; // boundaries
for j = 2:25; for k = 1:25; vm(j,k) = 0.5*(v(j,k)+v(j-1,k)); end; end;
for k = 1:25; vm(1,k) = vm(2,k);end; // boundaries
// eliminate grid points of low speeds
for j = 1:25; for k = 1:25;
speed(j,k) = sqrt(um(j,k)*um(j,k)+vm(j,k)*vm(j,k));
if speed(j,k) < 0.01; um(j,k) = 0.0; vm(j,k) = 0.0; end;
end; end;

drawlater; clf();

// 2d color plot of density field
Sgrayplot(x,y,rho,zminmax=[-0.1,0.01]);

// superimpose velocity arrows
champ(x,y,um,vm);

// specify graph & axis properties
a = gca();  a.font_size = 3; a.data_bounds = [0,0;50,50];
a.auto_ticks = ["off","off","on"]; a.sub_ticks = [3,3];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);

title("Time = "+string(0.01*int(100*time/24))+" days","fontsize",4,'position',[25 50]); // add title
xstring(22.5,2,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1;
xstring(2,22.5,"y (cm)");  // add z label
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
