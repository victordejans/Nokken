figure
[beta_vec,rho_min] = gen_fig_Kloomok_Muffley(0,20,1);
subplot(2,2,1)
plot(beta_vec, rho_min);
grid
xlabel('Beta (degrees)')
ylabel('Rho min (mm)')
title('Section 1')
[beta_vec,rho_min] = gen_fig_Kloomok_Muffley(20,35,2);
subplot(2,2,2)
plot(beta_vec, rho_min);
grid
xlabel('Beta (degrees)')
ylabel('Rho min (mm)')
title('Section 2')
[beta_vec,rho_min] = gen_fig_Kloomok_Muffley(0,0,1);
subplot(2,2,3)
plot(beta_vec, rho_min);
grid
xlabel('Beta (degrees)')
ylabel('Rho min (mm)')
title('Section 3 and 5 (dwells)')
[beta_vec,rho_min] = gen_fig_Kloomok_Muffley(35,0,6);
subplot(2,2,4)
plot(beta_vec, rho_min);
grid
xlabel('Beta (degrees)')
ylabel('Rho min (mm)')
title('Section 4')
