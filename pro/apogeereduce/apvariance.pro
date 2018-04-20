; variance returns the variance of data from Poisson + readout noise
;  gain and readout noise in electrons must be set

function apvariance,data,nread,gain,rn

if n_elements(gain) le 0 then gain=1.9
if n_elements(rn) le 0 then rn=18
rng=rn/gain
var=data/gain
neg=where(var le 0)
if neg[0] ge 0 then var[neg]=0
var+=nread*rng^2
return,var

end


