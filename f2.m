function y=f2(delta1,delta2,delta2d,beta_L,beta_C,phi_b,ib)
    y=-(delta2-delta1-2*pi*phi_b)/(pi*beta_L*beta_C)-sin(delta2)/beta_C-delta2d/beta_C+ib/beta_C;
end