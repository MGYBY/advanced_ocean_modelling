clear;

ans =messagebox("This calculator computes either phase speed ...
or period of surface gravity waves from user-specified values of total water depth ...
and wavelength. Do you want to continue?", "modal", "info", ["Yes" "No"]);

if ans == 1 then


k = x_choose_modeless(['Phase speed';'Wave period'],'What do you want to calculate?');
k
mode(-1); // silent mode
desc = "Enter values of total water depth (m) and wavelength (m)"
labls = ["Total water depth, H (m)"; "Wavelength, L (m)"];
typelst = list("vec",1,"vec",1);
deflts = ["1000";"100"];

//Call 'getvalue' to load values
[ok,H,L] = getvalue(desc,labls,typelst,deflts);

if ok then

exec("dispersion_relation.sce",-1) // load function

c = 0.0; t = 0.0;

if L*H > 0;
  [c,t] = disprel(L,H);
end;

if k == 1 then
Lines = ['For a total water depth of '+ string(H) + ' metres'; ... 
'and a wavelength of '+ string(L) + ' metres'; ...
'the phase speed of surface gravity waves is'; ...
'c = ' + string(int(100*c)/100) + ' metres per second'];
end;

if k == 2 then
Lines = ['For a total water depth of '+ string(H) + ' metres'; ...
'and a wavelength of '+ string(L) + ' metres,'; ...
'the period of surface gravity waves is'; ...
'T = ' + string(int(100*t)/100) + ' seconds'];
end;
messagebox(Lines);
else
  messagebox("No input available")
end;  

end;

