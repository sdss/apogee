function func_lincomb,x,y,par,im1=im1,im2=im2

model = im1*par[0] + im2*par[1] + par[2]

return,model
  
end
