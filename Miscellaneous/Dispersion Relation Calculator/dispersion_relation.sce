// this function calculates the phase speed c (m/s) and period t (s)
// of surface gravity waves
// for given wavelength L (m) and water depth H (m)
  function [c,t] = disprel(L,H)
    c =  sqrt(9.81*L/(2.*%pi)*tanh(2.*%pi*H/L));
    t = L/c;
  endfunction
