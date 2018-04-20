function func_poly2d,x,y,p

case n_elements(p) of
1:  a = p[0] 
3:  a = p[0] + p[1]*x + p[2]*y 
4:  a = p[0] + p[1]*x + p[2]*x*y + p[3]*y 
6:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x*y + p[4]*y + p[5]*y^2
8:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x*y + p[4]*(x^2)*y + $
        p[5]*x*y^2 + p[6]*y + p[7]*y^2
10:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x^3 + p[4]*x*y + p[5]*(x^2)*y + $
        p[6]*x*y^2 + p[7]*y + p[8]*y^2 + p[9]*y^3
11:  a = p[0] + p[1]*x + p[2]*x^2.0 + p[3]*x^3.0 + p[4]*x*y + p[5]*(x^2.0)*y + $
        p[6]*x*y^2.0 + p[7]*(x^2.0)*(y^2.0) + p[8]*y + p[9]*y^2.0 + p[10]*y^3.0
15:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x^3 + p[4]*x^4 + p[5]*y + p[6]*x*y + $
         p[7]*(x^2)*y + p[8]*(x^3)*y + p[9]*y^2 + p[10]*x*y^2 + p[11]*(x^2)*y^2 + $
         p[12]*y^3 + p[13]*x*y^3 + p[14]*y^4
21:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x^3 + p[4]*x^4 + p[5]*x^5 + p[6]*y + p[7]*x*y + $
         p[8]*(x^2)*y + p[9]*(x^3)*y + p[10]*(x^4)*y + p[11]*y^2 + p[12]*x*y^2 + $
         p[13]*(x^2)*y^2 + p[14]*(x^3)*y^2 + p[15]*y^3 + p[16]*x*y^3 + p[17]*(x^2)*y^3 + $ 
         p[18]*y^4 + p[19]*x*y^4 + p[20]*y^5
else: begin
       print,'Only 3, 4, 6, 8, 11 amd 15 parameters supported'
       return,-1
     end
endcase

return,a

end
