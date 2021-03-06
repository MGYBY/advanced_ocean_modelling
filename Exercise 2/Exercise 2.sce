//======================================
// Exercise 2: the bottom Ekman layer
//======================================
//Author: Jochen Kaempf, 2015 (update)

f = gcf(); scf(0); f.figure_size = [400,400];

in = read("uvprof1.dat",-1,3); // read data file
nz = 500; z = in(1:nz,1); u = in(1:nz,2); v = in(1:nz,3);
time = 0.0; dt = 2.0;
x = z; y = z; x(1:nz) = 0; y(1:nz) = 0;

for n = 1:150 // animation loop

time = n*dt;

// predict displacement of virtual floats
for i = 1:nz; x(i) = x(i)+dt*u(i); y(i) = y(i)+dt*v(i); end;

drawlater; clf; //Author: Jochen Kaempf, 2015 (update)

// draw graph frame
plot2d(0,0,-1,"031"," ",[-30,-20,30,40]);
ax = gca(); ax.font_size = 3;

// draw floats
nzz = 100; 
for i = 1:nz
  r = i/nz*1.0; r = max(r,0.3); d = r/2;
  xfarc(x(i)-d,y(i)+d,r,r,0,360*64) 
  p1=gce(); //get handle on current entity 
  p1.mark_foreground=0; p1.thickness=1; p1.fill_mode = "on"; p1.line_mode = "on";
end;

// draw flow vectors every 5 m
for i = 460:5:500
 xarrows([0 x(i)],[0 y(i)],40,1)
 p3=gce(); p3.mark_foreground=1; p3.thickness=1; p3.arrow_size=20;
 p3.segs_color = 1 ; 
end;

// draw displacement by ambient flow
xarrows([0 0],[0 y(1)],100,5); ss = min(time,100)/100;
p2=gce(); p2.mark_foreground=4; p2.thickness=3; p2.arrow_size=ss*40;

xgrid(1); // draw grid
p4=gce(); p4.mark_foreground=1; p4.thickness=1;

xstring(-5, -19,"x (m)");  // draw x label
txt=gce(); txt.font_size = 3;
xstring(-29, 0,"y (m)");  // draw y label
txt=gce(); txt.font_size = 3;

title("Time = "+string(int(time))+" secs","fontsize",3); // draw title

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

end; // end reference for animation loop

