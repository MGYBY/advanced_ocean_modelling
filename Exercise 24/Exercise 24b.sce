//=======================================
// Exercise 24: The abyssal circulation
//=======================================
// Animation of float movements
// Author: Jochen Kaempf, March 2015 (update)

f = gcf(); f.color_map = hsvcolormap(32); f.figure_size = [700,700]; scf(0);

// read input data
h0=read("h.dat",-1,51); 
trx1=read("TRx.dat",-1,5000); try1=read("TRy.dat",-1,5000); trz1=read("TRz.dat",-1,5000);
x = (0:20:1000)'; y = (0:20:1000)';  
[ntot ntra] =size(trx1);
time = 0.0;

// manipulate domain for graphical display
h1 = h0;
for j = 1:51; h1(j,51) = 1000; end;
for k = 1:51; h1(41,k) = 1000; end;
for j = 1:51; 
h1(j,31:32) = 1000;
h1(j,20:21) = 1000;
end;

mark(1:5000) = 0;
dummy = try1(1,1:ntra);
// mark floats according to different depth ranges
for i = 1:ntra
  if dummy(1,i) < 400 then 
      mark(i) = 1; 
  elseif dummy(1,i) > 800 then
      mark(i) = 2; 
  else
      mark(i) = 3; 
   end;
end;

//*******************************
for n = 1:ntot // animation loop
//*******************************

time = n-1; 

// grab data row
trax = trx1(n,1:ntra); tray = try1(n,1:ntra); traz = trz1(n,1:ntra);

drawlater; clf;

plot3d(x,y,-h1', 320,20," ",[0 3 4],[0,1000,0,1000,-1000,0]);

// draw floats

n1 = 0; n2 = 0; 
// sorting
for i = 1:ntra
  if mark(i) == 2 then
     n1 = n1+1; G1x(n1) = trax(1,i); G1y(n1) = tray(1,i); G1z(n1) = traz(1,i); 
  end;  
  if mark(i) == 1 then
    if traz(1,i) < 400 then // only display floats located in the upper 400 m
       n2 = n2+1; G2x(n2) = trax(1,i); G2y(n2) = tray(1,i); G2z(n2) = traz(1,i); 
    end;
  end;  
end;

//n2 = 0;

// draw first group
if n1 > 0 then
  c1(1:n1) = -9; 
  param3d1(G1x,G1y,list(-G1z,c1(1,:)),320, 20);
  f1 = gce(); f1.clip_state = "off"; f1.mark_size = 1; f1.mark_foreground = -1; f1.mark_background = 11; 
end;

// draw second group
if n2 > 0 then 
  c2(1:n2) = -14;
  param3d1(G2x,G2y,list(-G2z,c2(1,:)),320, 20);
  f2 = gce(); f2.clip_state = "off"; f2.mark_size = 2;  f2.mark_foreground = -1; f2.mark_background = 30; 
end;

clear G1x; clear G1y; clear G1z; clear c1;
clear G2x; clear G2y; clear G2z; clear c2;

// specify graph & axis properties
a = gca(); a.font_size = 3; a.data_bounds = [0,0,-1000;1000,1000,0];
a.auto_ticks = ["off","off","off"]; a.sub_ticks = [3,3,4];
a.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.y_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
a.z_ticks = tlist(["ticks", "locations","labels"],..
 [-1000 -500 0], ["-1000" "-500" "0"]);

xstring(730,-400,"x (m)");  // add x label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1; txt.clip_state = "off";
xstring(1400,210,"y (cm)");  // add y label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1; txt.clip_state = "off";
xstring(50,-250,"z (cm)");  // add z label
txt=gce(); txt.font_size = 4; txt.font_foreground = -1; txt.clip_state = "off";

title("Time = "+string(0.01*int(100*time))+" days","fontsize",4,'position',[350 1020]); // add title

drawnow;

// save frames as sequential GIF files (optional)
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
