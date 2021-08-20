//==========================================
// Exercise 3: Short Surface Gravity Waves
//==========================================
// Animation of equivalent vertical displacements of pressure surfaces
//Author: Jochen Kaempf, 2015 (update)

f = gcf(); scf(0); f.figure_size = [700,400]; f.children.font_size = 3;

// read input data
eta1=read("eta.dat",-1,101); dp1=read("dp.dat",-1,101); 
[ntot nx] = size(eta1); x = (0:5:500)'; 

for n = 1:100// animation loop

time = n; // time in seconds

//grab data blocks
itop = (n-1)*51+1; ibot = itop+50; 
dp = dp1(itop:ibot,1:101)'; eta = eta1(n,1:101)'; 

drawlater; clf();

// draw graphs
plot2d(x,5*eta,5); p1=gce(); p1.children.thickness=2;

for i = 1:26
  plot2d(x,5*dp(:,i)+1-i*2,2);//,'019','',[0 -40 500 10],[1,6,1,6]);
  p2=gce(); p2.children.thickness=1;
  b = gca(); b.font_size = 3; b.data_bounds = [0,-40;500,10];
  b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
  b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 100 200 300 400 500], ["0" "100" "200" "300" "400" "500"]);
  b.y_ticks = tlist(["ticks", "locations","labels"],..
  [-40 -30 -20 -10 0 10], ["-40" "-30" "-20" "-10" "0" "10"]);
    
end;

title("Time = "+string(int(time))+" secs","fontsize",3); // draw title
 
xstring(234, -38,"x (m)");  // draw x label
txt=gce(); txt.font_size = 4;
xstring(2, -22,"z (m)");  // draw z label
txt=gce(); txt.font_size = 4;

drawnow;

// save frames as sequential GIF files
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end // end reference for animation loop

