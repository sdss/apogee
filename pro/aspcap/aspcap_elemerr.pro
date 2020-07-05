function aspcap_elemerr,soln,te,fe,sn, quad=quad

  out=soln[0]+soln[1]*te+soln[2]*sn+soln[3]*fe
  sz=size(soln,/dim)
  if sz eq 5 then out+=soln[4]*te^2
  return,exp(out)

end

