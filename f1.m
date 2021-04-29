function y=f1(delta1,delta2,delta1d,beta_L,beta_C,phi_b)
    y=(delta2-delta1-2*pi*phi_b)/(pi*beta_L*beta_C)-sin(delta1)/beta_C-delta1d/beta_C;
end