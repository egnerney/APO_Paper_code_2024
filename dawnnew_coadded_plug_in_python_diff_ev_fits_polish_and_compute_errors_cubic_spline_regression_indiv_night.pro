
function given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now,Tec_a,nec_a,nsp_a,ns2p_a,nop_a,rho_max,rho_min, drho, yptsi_,nel_,tec_

  num_emiss = 8
  num_elements = n_elements(nec_a)
  epsilon_a = dblarr(num_elements,num_emiss)
  ;yptsi_,nel_,tec_,xwavi,yname

  ; stop
  triangulate,nel_,tec_,tr
  for j=0,num_emiss -1 do begin
    epsilon_a(*,j) = griddata(  nel_, tec_,reform(yptsi_(*,j)),xout = nec_a, yout = tec_a,/linear,triangles = tr )

  Endfor


  epsilon_a(*,0) = 1d-6*ns2p_a*epsilon_a(*,0)
  epsilon_a(*,1) = 1d-6*nop_a*epsilon_a(*,1)
  epsilon_a(*,2) = 1d-6*nop_a*epsilon_a(*,2)
  epsilon_a(*,3) = 1d-6*nsp_a*epsilon_a(*,3)
  epsilon_a(*,4) = 1d-6*nsp_a*epsilon_a(*,4)
  epsilon_a(*,5) = 1d-6*ns2p_a*epsilon_a(*,5)
  epsilon_a(*,6) = 1d-6*nsp_a*epsilon_a(*,6)
  epsilon_a(*,7) = 1d-6*nsp_a*epsilon_a(*,7)

  ; rho_c min of line of sight must be decreasing for this algorithm for increasing index fyi
  ;epsilon_a = dblarr(num_elements,8)

  rayleighs = dblarr(num_elements,8)

  ;sigma_n = dblarr(num_elements)


  conversion = 7.1492d9 ; 1RJ  in cm
  ; 4.6 - 7 finite tor
  n_elem_s = floor(2d*rho_max/drho) + 1
  s_LOS = drho*dindgen(n_elem_s)

  n_elem_x0_LOS = num_elements;floor((rho_max - rho_min)/drho) + 1
  y0 = - rho_max
  x_LOS = drho*dindgen(n_elem_x0_LOS) + rho_min ; x = x_0
  y_LOS = y0 + s_los

  rho_LOS = dblarr(n_elem_x0_LOS,n_elem_s);sqrt( y_LOS^2d + x_LOS^2d )
  idx_wantz = dblarr(n_elem_x0_LOS,n_elem_s)

  x_los = reverse(x_los)

  for i = 0 , n_elem_x0_LOS - 1 do rho_los(i,*) = sqrt( x_LOS(i)^2d + y_LOS^2d )

  for idx = 0,   num_elements - 1 do begin
    idx_want = where((reform(rho_LOS(idx,*)) le rho_max ) and (reform(rho_LOS(idx,*)) ge rho_min ))
    if ((idx_want(0) ge 0) and (n_elements( idx_want) ge 2 )) then begin
      ;where(x_los eq rho_LOS(idx,*))
      for linenum = 0, 7 do begin
        emiss_s = interpol(reform(epsilon_a(*,linenum)),x_LOS,reform(rho_LOS(idx,idx_want))) ; reform(epsilon_a(idx_want,linenum))
        rayleighs(idx,linenum) = conversion*INT_TABULATED( s_LOS(idx_want),emiss_s, /DOUBLE)
      endfor
    endif else if ((n_elements(idx_want) eq 1)  and (idx_want(0) ge 0)) then begin
    if (idx_want(0) + 1 lt n_elements(s_LOS)) then begin
        path_length = s_LOS[idx_want[0] + 1] - s_LOS[idx_want[0]]
    endif else begin
        path_length = drho
    endelse
    for linenum = 0, 7 do begin
      emiss_value = interpol(reform(epsilon_a(*,linenum)), x_LOS , reform(rho_LOS(idx,idx_want)))
      rayleighs(idx,linenum) = conversion * emiss_value * path_length
    endfor
    endif else begin
      rayleighs(idx,*) = replicate(0d,8)
    endelse


  endfor

  return, rayleighs
  ;sigma_n is given as optional output via keyword argument above

end

function given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now_float,Tec_a,nec_a,nsp_a,ns2p_a,nop_a,rho_max,rho_min, drho, yptsi_,nel_,tec_

  num_emiss = 8
  num_elements = n_elements(nec_a)
  epsilon_a = fltarr(num_elements,num_emiss)
  ;yptsi_,nel_,tec_,xwavi,yname

  ; stop
  for j=0,num_emiss -1 do begin
    triangulate,nel_,tec_,tr
    epsilon_a(*,j) = griddata(  nel_, tec_,reform(yptsi_(*,j)),xout = nec_a, yout = tec_a,/linear,triangles = tr )

  Endfor


  epsilon_a(*,0) = 1e-6*ns2p_a*epsilon_a(*,0)
  epsilon_a(*,1) = 1e-6*nop_a*epsilon_a(*,1)
  epsilon_a(*,2) = 1e-6*nop_a*epsilon_a(*,2)
  epsilon_a(*,3) = 1e-6*nsp_a*epsilon_a(*,3)
  epsilon_a(*,4) = 1e-6*nsp_a*epsilon_a(*,4)
  epsilon_a(*,5) = 1e-6*ns2p_a*epsilon_a(*,5)
  epsilon_a(*,6) = 1e-6*nsp_a*epsilon_a(*,6)
  epsilon_a(*,7) = 1e-6*nsp_a*epsilon_a(*,7)

  ; rho_c min of line of sight must be decreasing for this algorithm for increasing index fyi
  ;epsilon_a = dblarr(num_elements,8)

  rayleighs = fltarr(num_elements,8)

  ;sigma_n = dblarr(num_elements)


  conversion = 7.1492e9 ; 1RJ  in cm
  ; 4.6 - 7 finite tor
  n_elem_s = floor(2.*rho_max/drho) + 1
  s_LOS = drho*findgen(n_elem_s)

  n_elem_x0_LOS = floor((rho_max - rho_min)/drho) + 1
  y0 = - rho_max
  x_LOS = drho*findgen(n_elem_x0_LOS) + rho_min ; x = x_0
  y_LOS = y0 + s_los

  rho_LOS = fltarr(n_elem_x0_LOS,n_elem_s);sqrt( y_LOS^2d + x_LOS^2d )
  idx_wantz = fltarr(n_elem_x0_LOS,n_elem_s)

  x_los = reverse(x_los)

  for i = 0 , n_elem_x0_LOS - 1 do rho_los(i,*) = sqrt( x_LOS(i)^2. + y_LOS^2. )

  for idx = 0,   num_elements - 1 do begin
    idx_want = where((reform(rho_LOS(idx,*)) le rho_max ) and (reform(rho_LOS(idx,*)) ge rho_min ))
    if ((idx_want(0) ge 0) and (n_elements( idx_want) ge 2 )) then begin


      ;where(x_los eq rho_LOS(idx,*))
      for linenum = 0, 7 do begin
        emiss_s = interpol(reform(epsilon_a(*,linenum)),x_LOS,reform(rho_LOS(idx,idx_want))) ; reform(epsilon_a(idx_want,linenum))
        rayleighs(idx,linenum) = conversion*INT_TABULATED( s_LOS(idx_want),emiss_s)
      endfor
    endif else begin
      rayleighs(idx,*) = replicate(0.,8)
    endelse


  endfor

  return, rayleighs
  ;sigma_n is given as optional output via keyword argument above

end





function apo_deviates_3par_tec_pl_dusk,p,x=x,y=y,err=err
   COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_,rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  nec_funcform_test = dblarr(n_elements(x_grid))
  nsp_funcform_test = dblarr(n_elements(x_grid))
  ns2p_funcform_test = dblarr(n_elements(x_grid))
  nop_funcform_test = dblarr(n_elements(x_grid))

  ;nsp_to_nec_ratio = dblarr(n_elements(x_grid))

  ;flagsp = replicate(0d,n_elements(x_grid))
  
  flags2p = replicate(0d,n_elements(x_grid))
  ;flagop = replicate(0d,n_elements(x_grid))
 ; rhocne = p[0];5.27
  rhoRne =  rhoRne_ ;p[1];5.66d;p[1];5.69
  rhoORne = rhoORne_;5.82d;p[1];5.8 rhoRne_,rhoORne_,rhocne_
  ;A = p[1];1960
  ;Bin = p[2];0.2
  ;Bout = p[3];0.2
  ;C = p[4];3430
  ;Din = p[5];0.2
 ; Dout = p[6];0.2
 ; F = nec_pl;5.4d;p[14];5.4d

 ; E = A*Exp(-((rhoORne - rhocne)/Bout)^2d) + C*Exp(-((rhoORne - rhoRne)/Dout)^2d)
  rhocnsp = p[0] ;5.27
  rhoRnsp =  rhoRnsp_ ; 5.92d ;rhoRne - 0.09d;p[8];5.60
  rhoORnsp = rhoORne_;p[9] ;5.7
  G = p[1];1298
  Hin = p[2];0.2
  Hout = p[3];
  I = p[4];743
  Jin = p[5];0.08
  Jout = p[6];
  L = 3.37d + nec_pl

  K = G*Exp(-((rhoORnsp - rhocnsp)/Hout)^2d) + I*Exp(-((rhoORnsp - rhoRnsp)/Jout)^2d)
  
  
  ;rhocns2p = 5.05d
  rhoRns2p =  rhoRns2p_;5.7d
  rhoORns2p = rhoORne_;5.82d
  ;A = 1960d
  ;Bin = 0.3d
  ;Bout = 0.3d
  C2 = p[7];520d
  Din2 = p[8];0.2d
  Dout2 = p[9];0.2d
  F2 = nec_pl  + 1.13d

  E2 =  C2*Exp(-((rhoORns2p - rhoRns2p)/Dout2)^2d)




  rhocnop =  rhocnsp ;5.27d
  rhoRnop = rhoRne_;rhoRne - 0.09d;5.70d
  rhoORnop = rhoORne_;5.8d
  G2 = p[10];470d;1298d
  Hin2 = Hin ;0.2d
  Hout2 = Hout ;0.2d
  I2 = p[11];1106d;743d
  Jin2 = p[12];0.1d
  Jout2 = p[13];0.1d
  L2 = nec_pl + 0.96d

  K2 = G2*Exp(-((rhoORnop - rhocnop)/Hout2)^2d) + I2*Exp(-((rhoORnop - rhoRnop)/Jout2)^2d)
  
  tec = dblarr(n_elements(x_grid))
  A_TE1 = p[14]
  B_te1 = p[15]
  B_te2 = p[16]
  
  A_TE2 = p[17];2.1d;5.d; 3.07d


  beta_ = Alog(A_TE2 / A_TE1)/Alog(  rhoORne_ /   rho4069peak)



  for idx = 0, n_elements(x_grid) -1 do begin

    ;nec_funcform_test(idx) = E*(x_grid(idx)/rhoORne)^(-F)

    nsp_funcform_test(idx) = K*(x_grid(idx)/rhoORnsp)^(-L)

   ; if (x_grid(idx) lt rhoORne) then  nec_funcform_test(idx) = A*Exp(-((x_grid(idx) - rhocne)/Bout)^2d) + C*Exp(-((x_grid(idx) - rhoRne)/Dout)^2d)

    if (x_grid(idx) lt rhoORnsp) then  nsp_funcform_test(idx) = G*Exp(-((x_grid(idx) - rhocnsp)/Hout)^2d) + I*Exp(-((x_grid(idx) - rhoRnsp)/Jout)^2d)

    ;if (x_grid(idx) lt rhoRne) then  nec_funcform_test(idx) = A*Exp(-((x_grid(idx) - rhocne)/Bout)^2d) + C*Exp(-((x_grid(idx) - rhoRne)/Din)^2d)

    if (x_grid(idx) lt rhoRnsp) then  nsp_funcform_test(idx) = G*Exp(-((x_grid(idx) - rhocnsp)/Hout)^2d) + I*Exp(-((x_grid(idx) - rhoRnsp)/Jin)^2d)

   ; if (x_grid(idx) lt rhocne) then  nec_funcform_test(idx) = A*Exp(-((x_grid(idx) - rhocne)/Bin)^2d) + C*Exp(-((x_grid(idx) - rhoRne)/Din)^2d)

    if (x_grid(idx) lt rhocnsp) then  nsp_funcform_test(idx) = G*Exp(-((x_grid(idx) - rhocnsp)/Hin)^2d) + I*Exp(-((x_grid(idx) - rhoRnsp)/Jin)^2d)

    ;if ((nsp_funcform_test(idx) / nec_funcform_test(idx)) gt 0.8d) then flagsp(idx) = 1d

    ns2p_funcform_test(idx) = E2*(x_grid(idx)/rhoORns2p)^(-F2)

    nop_funcform_test(idx) = K2*(x_grid(idx)/rhoORnop)^(-L2)

    if (x_grid(idx) lt rhoORns2p) then  ns2p_funcform_test(idx) = C2*Exp(-((x_grid(idx) - rhoRns2p)/Dout2)^2d)

    if (x_grid(idx) lt rhoORnop) then  nop_funcform_test(idx) = G2*Exp(-((x_grid(idx) - rhocnop)/Hout2)^2d) + I2*Exp(-((x_grid(idx) - rhoRnop)/Jout2)^2d)

    if (x_grid(idx) lt rhoRns2p) then  ns2p_funcform_test(idx) =  C2*Exp(-((x_grid(idx) - rhoRns2p)/Din2)^2d)

    if (x_grid(idx) lt rhoRnop) then  nop_funcform_test(idx) = G2*Exp(-((x_grid(idx) - rhocnop)/Hout2)^2d) + I2*Exp(-((x_grid(idx) - rhoRnop)/Jin2)^2d)

    if (x_grid(idx) lt rhocnop) then  nop_funcform_test(idx) = G2*Exp(-((x_grid(idx) - rhocnop)/Hin2)^2d) + I2*Exp(-((x_grid(idx) - rhoRnop)/Jin2)^2d)
   ; if ((ns2p_funcform_test(idx) / nec_funcform_test(idx)) gt 0.4d) then flags2p(idx) = 1d
   ; if ((nop_funcform_test(idx) / nec_funcform_test(idx)) gt 0.8d) then flagop(idx) = 1d
   nec_funcform_test(idx) = (nsp_funcform_test(idx)  + 2d*ns2p_funcform_test(idx)  + nop_funcform_test(idx) )/(totchmix_for_apo(idx))

if (x_grid(idx) le  rho4069peak ) then tec[idx] = A_Te1*((x_grid[idx] / rho4069peak  )^(B_te1))
if (x_grid(idx) gt  rho4069peak ) then tec[idx] = A_Te1*((x_grid[idx] / rho4069peak  )^(beta_))
if (x_grid(idx) gt rhoORne_ ) then tec[idx] = A_Te2*((x_grid[idx] / rhoORne_ )^(B_te2))

   
   if (((ns2p_funcform_test(idx) / nec_funcform_test(idx)) gt 0.05d) and (x_grid(idx) le 5.3d)) then flags2p(idx) = 1d

   if ( nec_funcform_test(idx)  lt 0.000001d ) then begin
     nsp_mix_temp = nsp_funcform_test(idx)/nec_funcform_test(idx)
     ns2p_mix_temp = ns2p_funcform_test(idx)/nec_funcform_test(idx)
     nop_mix_temp = nop_funcform_test(idx)/nec_funcform_test(idx)
     nec_funcform_test(idx)  = 0.000001d
     nsp_funcform_test(idx)  = nsp_mix_temp*nec_funcform_test(idx)
     ns2p_funcform_test(idx)  = ns2p_mix_temp*nec_funcform_test(idx)
     nop_funcform_test(idx)  = nop_mix_temp*nec_funcform_test(idx)
   endif



  endfor
  
  
 




 ; flagtotsp = total(flagsp,/double)
  
  flagtots2p = total(flags2p,/double)
  ;flagtotop = total(flagop,/double)

  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(tec,nec_funcform_test,nsp_funcform_test,ns2p_funcform_test,nop_funcform_test,max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)


  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
  idx_want =  [0,1,2,3,5,6,7] ;
  model = model[*,idx_want]

  model = reform(model,n_elements(model))
  
  penalt=1d6
   buffer_val = 0d
   ;model = model
   
 ; if (flagtotsp gt 0) then model = model + penalt*flagtotsp/(double(n_elements(model))) ; penalty for neutrality condition broken 
  ;if ( rhoRne  ge  rhoORne) then model = model + penalt/(double(n_elements(model)))
  ;if ( rhoRnsp ge  rhoORnsp) then model = model + penalt/(double(n_elements(model)))
  ;if ( rhoRnsp ge  rhoRne) then model = model + penalt/(double(n_elements(model)))
  if (flagtots2p gt 0) then model = model + penalt*flagtots2p/(double(n_elements(model))) ; penalty for neutrality condition broken
  ;if (flagtotop gt 0) then model = model + penalt*flagtotop/(double(n_elements(model)))
  
  result = (y-model)/err

  ; if (G > A) then model = model + 1d6*(1d - (A - G)) ; penalty for no physical cold torus
  ;y needs to be setup to for [model] =  reform(dblarr(num_elements,8),num_elements*8) ; well 6 now
  RETURN, result

end





function apo_deviates_emission_table_allside_simult_foursix_to_seven_7lines_no_func_form_fixed_necbyratio_tec_pl_float_dawn,p,x=x,y=y,err=err

  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_, rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  ; nec_tempi  = p[0:nx_grid - 1]

  nsp_tempi = p[0:nx_grid - 1]

  ns2p_tempi = p[nx_grid:2*nx_grid - 1]

  nop_tempi = p[2*nx_grid:3*nx_grid - 1]

  tec_tempi = p[3*nx_grid:4*nx_grid - 1]

  nec_tempi  =  (nsp_tempi  + 2.*ns2p_tempi  + nop_tempi )/(totchmix_for_apo)
  ;nsp_to_nec_tot = total(nsp_tempi,/double)/total(nec_tempi,/double)
  ;ns2p_to_nec_tot = total(ns2p_tempi,/double)/total(nec_tempi,/double)
  ; nop_to_nec_tot = total(nop_tempi,/double)/total(nec_tempi,/double)
  ;totchmix = total(nsp_tempi + 2d*nsp_tempi + nop_tempi)/total(nec_tempi)
  ; flaggsp = replicate(0d,nx_grid)
  ; flaggs2p = replicate(0d,nx_grid)
  ;flaggop = replicate(0d,nx_grid)
  ;for idxx = 0, nx_grid - 1 do begin
  ;   if ((nsp_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggsp(idxx) = 1d
  ; if (((ns2p_tempi(idxx)/nec_tempi(idxx)) ge 0.05d) and (x_grid(idxx) le 5.0d))  then flaggs2p(idxx) = 1d
  ;  if ((nop_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggop(idxx) = 1d
  ;   if ((x_grid(idxx) lt (rhoRne_ - 0.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.1d)) then flaggs2p(idxx) = 1d
  ;  if ((x_grid(idxx) lt (rhoRne_ - 1.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.01d)) then flaggs2p(idxx) = 1d
  ; endfor
  ; flaggtotsp = total(flaggsp,/double)
  ;flaggtots2p = total(flaggs2p,/double)
  ;flaggtotop = total(flaggop,/double)

  ;model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed(te_fixed,p[0:nx_grid - 1],p[nx_grid:2*nx_grid - 1],p[2*nx_grid:3*nx_grid - 1],p[3*nx_grid:4*nx_grid - 1],max(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now_float( tec_tempi,nec_tempi,nsp_tempi,ns2p_tempi,nop_tempi,float(max(x_grid)),float(min(x_grid)), float(x_step_size),float(yptsi_in_2dto1d),float(nel2dto1d),float(tec2dto1d))

  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
  idx_want = [0,1,2,3,5,6,7];[3,6,7];[1,2,3,5,6,7]
  model = model[*,idx_want]

  model = reform(model,n_elements(model))

  ;if (flaggtot gt 0) then model = model + 1d4*flaggtot/(double(n_elements(model)))
  ;model = float(model)
  ;penaltt=1e6
  ;buffer_val2 = 0d
  ;model = model
  ;if (flaggtotsp gt 0) then model = model + penaltt*flaggtotsp/(double(n_elements(model))) ; penalty for neutrality condition broken
  ;  if (flaggtots2p gt 0) then model = model + penaltt*flaggtots2p/(double(n_elements(model))) ; penalty for neutrality condition broken
  ; if (flaggtotop gt 0) then model = model + penaltt*flaggtotop/(double(n_elements(model)))

  ; help,y
  ;help,err
  ;help,model
  ;stop
  result = (y - model)/err

  RETURN, result

end

function apo_deviates_emission_table_allside_simult_foursix_to_seven_7lines_no_func_form_fixed_necbyratio_tec_pl_double_dawn,p,x=x,y=y,err=err

  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_ , rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  ; nec_tempi  = p[0:nx_grid - 1]

  nsp_tempi = p[0:nx_grid - 1]

  ns2p_tempi = p[nx_grid:2*nx_grid - 1]

  nop_tempi = p[2*nx_grid:3*nx_grid - 1]

  tec_tempi = p[3*nx_grid:4*nx_grid - 1]

  nec_tempi  =  (nsp_tempi  + 2.d*ns2p_tempi  + nop_tempi )/ totchmix_for_apo 
  ;nsp_to_nec_tot = total(nsp_tempi,/double)/total(nec_tempi,/double)
  ;ns2p_to_nec_tot = total(ns2p_tempi,/double)/total(nec_tempi,/double)
  ; nop_to_nec_tot = total(nop_tempi,/double)/total(nec_tempi,/double)
  ;totchmix = total(nsp_tempi + 2d*nsp_tempi + nop_tempi)/total(nec_tempi)
  ; flaggsp = replicate(0d,nx_grid)
  ; flaggs2p = replicate(0d,nx_grid)
  ;flaggop = replicate(0d,nx_grid)
  ;for idxx = 0, nx_grid - 1 do begin
  ;   if ((nsp_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggsp(idxx) = 1d
  ; if (((ns2p_tempi(idxx)/nec_tempi(idxx)) ge 0.05d) and (x_grid(idxx) le 5.0d))  then flaggs2p(idxx) = 1d
  ;  if ((nop_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggop(idxx) = 1d
  ;   if ((x_grid(idxx) lt (rhoRne_ - 0.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.1d)) then flaggs2p(idxx) = 1d
  ;  if ((x_grid(idxx) lt (rhoRne_ - 1.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.01d)) then flaggs2p(idxx) = 1d
  ; endfor
  ; flaggtotsp = total(flaggsp,/double)
  ;flaggtots2p = total(flaggs2p,/double)
  ;flaggtotop = total(flaggop,/double)

  ;model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed(te_fixed,p[0:nx_grid - 1],p[nx_grid:2*nx_grid - 1],p[2*nx_grid:3*nx_grid - 1],p[3*nx_grid:4*nx_grid - 1],max(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now( tec_tempi,nec_tempi,nsp_tempi,ns2p_tempi,nop_tempi,max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
  idx_want = [0,1,2,3,5,6,7];[3,6,7];[1,2,3,5,6,7]
  model = model[*,idx_want]

  model = reform(model,n_elements(model))

  ;if (flaggtot gt 0) then model = model + 1d4*flaggtot/(double(n_elements(model)))
  ;model = float(model)
  ;penaltt=1e6
  ;buffer_val2 = 0d
  ;model = model
  ;if (flaggtotsp gt 0) then model = model + penaltt*flaggtotsp/(double(n_elements(model))) ; penalty for neutrality condition broken
  ;  if (flaggtots2p gt 0) then model = model + penaltt*flaggtots2p/(double(n_elements(model))) ; penalty for neutrality condition broken
  ; if (flaggtotop gt 0) then model = model + penaltt*flaggtotop/(double(n_elements(model)))

  ; help,y
  ;help,err
  ;help,model
  ;stop
  result = (y - model)/err

  RETURN, result

end

function apo_deviates_CSR_fixed_mixr,p,x=x,y=y,err=err

  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_ , rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  nec_tempi  = p[0:nx_grid - 1]

  nsp_tempi =  nec_tempi*nspmix

  ns2p_tempi = nec_tempi*ns2pmix

  nop_tempi = nec_tempi*nopmix

  tec_tempi = p[nx_grid:2*nx_grid - 1]
  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now( tec_tempi,nec_tempi,nsp_tempi,ns2p_tempi,nop_tempi,max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)


  idx_want = [0,1,2,3,5,6,7];[3,6,7];[1,2,3,5,6,7]
  model = model[*,idx_want]

  model = reform(model,n_elements(model))
  secondiv_nec = dblarr(nx_grid )
  secondiv_tec = secondiv_nec
  for i = 0,nx_grid - 1 do begin

if (i eq 0) then begin
    secondiv_nec[i] = (nec_tempi[i + 2] - 2d*nec_tempi[i + 1] + nec_tempi[i]) / x_step_size^2d
    secondiv_tec[i] = (tec_tempi[i + 2] - 2d*tec_tempi[i + 1] + tec_tempi[i]) / x_step_size^2d
end else if (i eq (nx_grid - 1)) then begin
    secondiv_nec[i] = (nec_tempi[i] - 2d*nec_tempi[i - 1] + nec_tempi[i - 2]) / x_step_size^2d
    secondiv_tec[i] = (tec_tempi[i] - 2d*tec_tempi[i - 1] + tec_tempi[i - 2]) / x_step_size^2d
end else begin
    secondiv_nec[i] = (nec_tempi[i + 1] - 2d*nec_tempi[i] + nec_tempi[i - 1]) / x_step_size^2d
    secondiv_tec[i] = (tec_tempi[i + 1] - 2d*tec_tempi[i] + tec_tempi[i - 1]) / x_step_size^2d
end
    
  endfor

  ;pen_1 = sqrt(lambdaa1*int_tabulated(x_grid,secondiv_nec^2d))
  
  ;pen_2 = sqrt(lambdaa2*int_tabulated(x_grid,secondiv_tec^2d))
  
  pen_1 = sqrt(lambdaa1*x_step_size)*(secondiv_nec)

  pen_2 = sqrt(lambdaa2*x_step_size)*(secondiv_tec)
  
  result = [(y[0:n_elements(model)-1] - model)/err[0:n_elements(model)-1], pen_1,pen_2]
  ;print,pen_1,pen_2
  ;stop
  RETURN, result

end

function apo_deviates_emission_table_allside_simult_foursix_to_seven_7lines_no_func_form_fixed_mixr_ratios_fit_nec_tec_pl_double_only,p,x=x,y=y,err=err

  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_, rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  nec_tempi  = p[0:nx_grid - 1]

  nsp_tempi =  nec_tempi*nspmix

  ns2p_tempi = nec_tempi*ns2pmix

  nop_tempi = nec_tempi*nopmix

  tec_tempi = p[nx_grid:2*nx_grid - 1]

  ;(nsp_tempi  + 2d*ns2p_tempi  + nop_tempi )/(totchmix_for_apo)
  ;nsp_to_nec_tot = total(nsp_tempi,/double)/total(nec_tempi,/double)
  ;ns2p_to_nec_tot = total(ns2p_tempi,/double)/total(nec_tempi,/double)
  ; nop_to_nec_tot = total(nop_tempi,/double)/total(nec_tempi,/double)
  ;totchmix = total(nsp_tempi + 2d*nsp_tempi + nop_tempi)/total(nec_tempi)
  ; flaggsp = replicate(0d,nx_grid)
  ; flaggs2p = replicate(0d,nx_grid)
  ;flaggop = replicate(0d,nx_grid)
  ;for idxx = 0, nx_grid - 1 do begin
  ;   if ((nsp_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggsp(idxx) = 1d
  ; if (((ns2p_tempi(idxx)/nec_tempi(idxx)) ge 0.05d) and (x_grid(idxx) le 5.0d))  then flaggs2p(idxx) = 1d
  ;  if ((nop_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggop(idxx) = 1d
  ;   if ((x_grid(idxx) lt (rhoRne_ - 0.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.1d)) then flaggs2p(idxx) = 1d
  ;  if ((x_grid(idxx) lt (rhoRne_ - 1.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.01d)) then flaggs2p(idxx) = 1d
  ; endfor
  ; flaggtotsp = total(flaggsp,/double)
  ;flaggtots2p = total(flaggs2p,/double)
  ;flaggtotop = total(flaggop,/double)

  ;model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed(te_fixed,p[0:nx_grid - 1],p[nx_grid:2*nx_grid - 1],p[2*nx_grid:3*nx_grid - 1],p[3*nx_grid:4*nx_grid - 1],max(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now( tec_tempi,nec_tempi,nsp_tempi,ns2p_tempi,nop_tempi,max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
  idx_want = [0,1,2,3,5,6,7];[3,6,7];[1,2,3,5,6,7]
  model = model[*,idx_want]

  model = reform(model,n_elements(model))

  ;if (flaggtot gt 0) then model = model + 1d4*flaggtot/(double(n_elements(model)))
  ;model = float(model)
  ;penaltt=1e6
  ;buffer_val2 = 0d
  ;model = model
  ;if (flaggtotsp gt 0) then model = model + penaltt*flaggtotsp/(double(n_elements(model))) ; penalty for neutrality condition broken
  ;  if (flaggtots2p gt 0) then model = model + penaltt*flaggtots2p/(double(n_elements(model))) ; penalty for neutrality condition broken
  ; if (flaggtotop gt 0) then model = model + penaltt*flaggtotop/(double(n_elements(model)))

  ; help,y
  ;help,err
  ;help,model
  ;stop
  result = (y - model)/err

  ;print,result

  ;stop

  RETURN, result

end

function apo_deviates_emission_table_allside_simult_foursix_to_seven_7lines_no_func_form_fixed_mixr_ratios_fit_nec_tec_pl_double_L2_regularization_dawn_only,p,x=x,y=y,err=err

  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_,rhoRns2p_ , rho4069peak, lambdaa1, lambdaa2
  nx_grid = n_elements(x_grid)
  nec_tempi  = p[0:nx_grid - 1]

  nsp_tempi =  nec_tempi*nspmix

  ns2p_tempi = nec_tempi*ns2pmix

  nop_tempi = nec_tempi*nopmix

  tec_tempi = p[nx_grid:2*nx_grid - 1]

  ;(nsp_tempi  + 2d*ns2p_tempi  + nop_tempi )/(totchmix_for_apo)
  ;nsp_to_nec_tot = total(nsp_tempi,/double)/total(nec_tempi,/double)
  ;ns2p_to_nec_tot = total(ns2p_tempi,/double)/total(nec_tempi,/double)
  ; nop_to_nec_tot = total(nop_tempi,/double)/total(nec_tempi,/double)
  ;totchmix = total(nsp_tempi + 2d*nsp_tempi + nop_tempi)/total(nec_tempi)
  ; flaggsp = replicate(0d,nx_grid)
  ; flaggs2p = replicate(0d,nx_grid)
  ;flaggop = replicate(0d,nx_grid)
  ;for idxx = 0, nx_grid - 1 do begin
  ;   if ((nsp_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggsp(idxx) = 1d
  ; if (((ns2p_tempi(idxx)/nec_tempi(idxx)) ge 0.05d) and (x_grid(idxx) le 5.0d))  then flaggs2p(idxx) = 1d
  ;  if ((nop_tempi(idxx)/nec_tempi(idxx)) gt 0.8d) then flaggop(idxx) = 1d
  ;   if ((x_grid(idxx) lt (rhoRne_ - 0.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.1d)) then flaggs2p(idxx) = 1d
  ;  if ((x_grid(idxx) lt (rhoRne_ - 1.5d)) and ((ns2p_tempi(idxx)/nec_tempi(idxx)) gt 0.01d)) then flaggs2p(idxx) = 1d
  ; endfor
  ; flaggtotsp = total(flaggsp,/double)
  ;flaggtots2p = total(flaggs2p,/double)
  ;flaggtotop = total(flaggop,/double)

  ;model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed(te_fixed,p[0:nx_grid - 1],p[nx_grid:2*nx_grid - 1],p[2*nx_grid:3*nx_grid - 1],p[3*nx_grid:4*nx_grid - 1],max(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
  model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now( tec_tempi,nec_tempi,nsp_tempi,ns2p_tempi,nop_tempi,max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
  idx_want = [0,1,2,3,5,6,7];[3,6,7];[1,2,3,5,6,7]
  model = model[*,idx_want]

  model = reform(model,n_elements(model))

  ;if (flaggtot gt 0) then model = model + 1d4*flaggtot/(double(n_elements(model)))
  ;model = float(model)
  ;penaltt=1e6
  ;buffer_val2 = 0d
  ;model = model
  ;if (flaggtotsp gt 0) then model = model + penaltt*flaggtotsp/(double(n_elements(model))) ; penalty for neutrality condition broken
  ;  if (flaggtots2p gt 0) then model = model + penaltt*flaggtots2p/(double(n_elements(model))) ; penalty for neutrality condition broken
  ; if (flaggtotop gt 0) then model = model + penaltt*flaggtotop/(double(n_elements(model)))

  ; help,y
  ;help,err
  ;help,model
  ;stop
  result1 = (y[0:n_elements(model)-1] - model)/err[0:n_elements(model)-1]

lambdum = total(result1^2d,/double)/10d
  ;lambdaa = 3.4d-7  ;8.5d-6;1.7d-6;0.000000017d;
  result = [result1,total(sqrt(lambdum)*p,/double)]

  ;print,result

  ;stop

  RETURN, result

end


function tanh_deviates,p,x=x,y=y,err=err


  model = -p[0]*tanh((x - p[1])/p[2]) + p[3]

  
  result = (y - model)/err

  RETURN, result

end


pro Dawnnew_coadded_plug_in_python_diff_ev_fits_polish_and_compute_errors_cubic_spline_regression_indiv_night

  restore, filename ='grids.sav', /verbose

  ;COMMON block1, i, n_in, Te_in, xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, x_grid, x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,te_noww,pp
  ;COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,te_fixed,nec_fixed_found
 ; COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed
 ; COMMON block2,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne,rhoORne,rhocne
  COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa,rhoRnsp_, rhoRns2p_ , rho4069peak, lambdaa1, lambdaa2 ;print,min(x_grid),max(x_grid)
  ;help,x_grid
  ;stop
  ;help,x_grid
  ;help,DAWN_GRID
  ;help,err_DAWN_GRIDk
  ;help,dusk_grid
  ;help,err_dusk_grid
  ;x_grid=reverse(x_grid)

  p1=errorplot(x_grid,dusk_grid(*,2),ERR_dusk_grid(*,2),NAME='OII (O+) 3729 Angstroms',ERRORBAR_COLOR="blue",color="blue",title='Nominal Coadded APO Dusk Profile From Carl',xtitle='Dusk Distance (RJ)',ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400])

  p2=errorplot(x_grid,dusk_grid(*,6),ERR_dusk_grid(*,6),NAME='SII (S+) 6716 Angstroms',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p3=errorplot(x_grid,dusk_grid(*,7),ERR_dusk_grid(*,7),NAME='SII (S+) 6731 Angstroms',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p4=errorplot(x_grid,dusk_grid(*,0),ERR_dusk_grid(*,0),NAME='SIII (S++) 3722 Angstroms',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p5=errorplot(x_grid,dusk_grid(*,1),ERR_dusk_grid(*,1),NAME='OII (O+) 3726 Angstroms',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p6=errorplot(x_grid,dusk_grid(*,5),ERR_dusk_grid(*,5),NAME='SIII (S++) 6312 Angstroms',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p7=errorplot(x_grid,dusk_grid(*,3),ERR_dusk_grid(*,3),NAME='SII (S+) 4069 Angstroms',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p8=errorplot(x_grid,dusk_grid(*,4),ERR_dusk_grid(*,4),NAME='SII (S+) 4076 Angstroms',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


  leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.9,300.], $
    /DATA, /AUTO_TEXT_COLOR)

  ;stop

  p1=errorplot(x_grid,dawn_grid(*,2),ERR_dawn_grid(*,2),NAME='OII (O+) 3729 Angstroms',ERRORBAR_COLOR="blue",color="blue",title='Nominal Coadded APO Dawn Profile From Carl',xtitle='Dawn Distance (RJ)',ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400],xrange=[7.5,4.5])

  p2=errorplot(x_grid,dawn_grid(*,6),ERR_dawn_grid(*,6),NAME='SII (S+) 6716 Angstroms',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p3=errorplot(x_grid,dawn_grid(*,7),ERR_dawn_grid(*,7),NAME='SII (S+) 6731 Angstroms',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p4=errorplot(x_grid,dawn_grid(*,0),ERR_dawn_grid(*,0),NAME='SIII (S++) 3722 Angstroms',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p5=errorplot(x_grid,dawn_grid(*,1),ERR_dawn_grid(*,1),NAME='OII (O+) 3726 Angstroms',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p6=errorplot(x_grid,dawn_grid(*,5),ERR_dawn_grid(*,5),NAME='SIII (S++) 6312 Angstroms',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p7=errorplot(x_grid,dawn_grid(*,3),ERR_dawn_grid(*,3),NAME='SII (S+) 4069 Angstroms',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p8=errorplot(x_grid,dawn_grid(*,4),ERR_dawn_grid(*,4),NAME='SII (S+) 4076 Angstroms',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


  leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.9,300.], $
    /DATA, /AUTO_TEXT_COLOR)

  ;stop

  ;;;;;;;;;;;;;;


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvalues=5

  xmin_want = 4.6d;5.
  xmax_want = 7.0d

  dx_want = 0.04d


  ;ninterpsteps = round((xinterpmax - xinterpmin)/xinterpstep) + 1
  factor = round(dx_want/0.02d) ; dx_want/og_dx_had
  n_want = round(n_elements(x_grid)/factor)
  idx_want = factor*indgen(n_want)
  ;x_want = ;xinterpstep*dindgen(ninterpsteps) + xinterpmin

  ;ninterpsteps = n_elements(x_want)


  ; interped_dawn_grid=dblarr(ninterpsteps,8)
  ; interped_dusk_grid=dblarr(ninterpsteps,8)
  ; interped_err_dawn_grid=dblarr(ninterpsteps,8)
  ; interped_err_dusk_grid=dblarr(ninterpsteps,8)


  dawn_grid_new = dblarr(n_elements(x_grid),8)
  dusk_grid_new = dblarr(n_elements(x_grid),8)

  err_dawn_grid_new = dblarr(n_elements(x_grid),8)
  err_dusk_grid_new = dblarr(n_elements(x_grid),8)

  for i=0,7 do begin
    for j = 2,n_elements(x_grid)-3 do begin

      ;5 points running average and associated errors/uncertainties assuming Gaussian

      dawn_grid_new[j,i] = (dawn_grid[j - 2,i] + dawn_grid[j - 1,i] + dawn_grid[j,i] + dawn_grid[j + 1,i] + dawn_grid[j + 2,i])/double(nvalues) ; ts_smooth(dawn_grid[*,i],nvalues)

      dusk_grid_new[j,i] =  (dusk_grid[j - 2,i] + dusk_grid[j - 1,i] + dusk_grid[j,i] + dusk_grid[j + 1,i] + dusk_grid[j + 2,i])/double(nvalues) ;ts_smooth(dusk_grid[*,i],nvalues)


      if (dawn_grid_new[j,i] lt 0d) then dawn_grid_new[j,i] = 0.1d

      if (dusk_grid_new[j,i] lt 0d) then dusk_grid_new[j,i] = 0.1d

      ;interped_dawn_grid[j,i] = ;interpol(reform(dawn_grid[*,i]),x_grid, x_interps)

      ;interped_dusk_grid[j,i] = ;interpol(dusk_grid[*,i],x_grid, x_interps)


      err_dawn_grid_new[j,i] = Sqrt( err_dawn_grid[j - 2,i]^2d + err_dawn_grid[j - 1,i]^2d + err_dawn_grid[j,i]^2d + err_dawn_grid[j + 1,i]^2d + err_dawn_grid[j + 2,i]^2d )/double(nvalues);ts_smooth(err_dawn_grid[*,i],nvalues)

      err_dusk_grid_new[j,i] =  Sqrt( err_dusk_grid[j - 2,i]^2d + err_dusk_grid[j - 1,i]^2d + err_dusk_grid[j,i]^2d + err_dusk_grid[j + 1,i]^2d + err_dusk_grid[j + 2,i]^2d )/double(nvalues);ts_smooth(err_dusk_grid[*,i],nvalues)



     ; if (dawn_grid_new[j,i] le 0d) then err_dawn_grid_new[j,i] = 10000d ; big number

      ;if (dusk_grid_new[j,i] le 0d) then err_dusk_grid_new[j,i] = 10000d


      ; interped_err_dawn_grid[j,i] = interpol(err_dawn_grid[*,i],x_grid, x_interps)

      ; interped_err_dusk_grid[j,i] = interpol(err_dusk_grid[*,i],x_grid, x_interps)
    endfor
  endfor

  dawn_grid = dawn_grid_new[idx_want,*] ;interped_dawn_grid
  dusk_grid = dusk_grid_new[idx_want,*];interped_dusk_grid

  err_dawn_grid = err_dawn_grid_new[idx_want,*];interped_err_dawn_grid/sqrt(double(nvalues))  ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg
  err_dusk_grid = err_dusk_grid_new[idx_want,*];interped_err_dusk_grid/sqrt(double(nvalues)) ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg



  x_grid = x_grid[idx_want];x_interps


  idx_want = where((x_grid ge xmin_want) and (x_grid le xmax_want))

  dawn_grid = dawn_grid[idx_want,*] ;interped_dawn_grid
  dusk_grid = dusk_grid[idx_want,*];interped_dusk_grid

  err_dawn_grid = err_dawn_grid[idx_want,*];interped_err_dawn_grid/sqrt(double(nvalues))  ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg
  err_dusk_grid = err_dusk_grid[idx_want,*];interped_err_dusk_grid/sqrt(double(nvalues)) ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg



  x_grid = x_grid[idx_want];x_interps




  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


 ; p1=errorplot(x_grid,dawn_grid(*,6)/dawn_grid(*,7),(dawn_grid(*,6)/dawn_grid(*,7))*sqrt((ERR_dawn_grid(*,6)/dawn_grid(*,6))^2d + (ERR_dawn_grid(*,7)/dawn_grid(*,7))^2d),NAME='6716Å to 6731Å Ratio',ERRORBAR_COLOR="indian_red",color="indian_red",title='Dawn 5 point Avg 0.04 Resolution',xtitle='Dawn $\rho_c$ ($R_J$)',ytitle='Line Ratio',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,1.5],layout=[1,2,1])

  ;  p2=errorplot(x_grid,dawn_grid(*,3)/dawn_grid(*,7),(dawn_grid(*,3)/dawn_grid(*,7))*sqrt((ERR_dawn_grid(*,3)/dawn_grid(*,3))^2d + (ERR_dawn_grid(*,7)/dawn_grid(*,7))^2d),NAME='4069Å to 6731Å Ratio',ERRORBAR_COLOR="peru",color="peru",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  ;p3=errorplot(x_grid,dusk_grid(*,6)/dusk_grid(*,7),(dusk_grid(*,6)/dusk_grid(*,7))*sqrt((ERR_dusk_grid(*,6)/dusk_grid(*,6))^2d + (ERR_dusk_grid(*,7)/dusk_grid(*,7))^2d),NAME='6716Å to 6731Å Ratio',ERRORBAR_COLOR="indian_red",color="indian_red",title='Dusk 5 point Avg 0.04 Resolution',xtitle='Dusk $\rho_c$ ($R_J$)',ytitle='Line Ratio',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,1.5],layout=[1,2,2],/current)

  ;  p4=errorplot(x_grid,dusk_grid(*,3)/dusk_grid(*,7),(dusk_grid(*,3)/dusk_grid(*,7))*sqrt((ERR_dusk_grid(*,3)/dusk_grid(*,3))^2d + (ERR_dusk_grid(*,7)/dusk_grid(*,7))^2d),NAME='4069Å to 6731Å Ratio',ERRORBAR_COLOR="peru",color="peru",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


 ; legday = LEGEND(TARGET=[p3,p4], POSITION=[7.1,2.4], $
  ;  /DATA, /AUTO_TEXT_COLOR)
    
   ; p3=errorplot(x_grid,dusk_grid(*,6)/dusk_grid(*,7),(dusk_grid(*,6)/dusk_grid(*,7))*sqrt((ERR_dusk_grid(*,6)/dusk_grid(*,6))^2d + (ERR_dusk_grid(*,7)/dusk_grid(*,7))^2d),NAME='6716Å to 6731Å Ratio',ERRORBAR_COLOR="green",color="green",title='Dusk 5 point Avg 0.04 Resolution',xtitle='Dusk $\rho_c$ ($R_J$)',ytitle='Line Ratio',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0.7,0.83],xrange=[5.4,5.85])


      ;p33=errorplot(x_grid,dusk_grid(*,7)*(0.82d/Max(dusk_grid(*,7))),err_dusk_grid(*,7)*(0.82d/Max(dusk_grid(*,7))),NAME='Normalized 6731Å',ERRORBAR_COLOR="brown",color="brown",title='Dusk 5 point Avg 0.04 Resolution',xtitle='Dusk $\rho_c$ ($R_J$)',ytitle='',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,xrange=[5.4,5.85],yrange=[0.7,0.83],/overplot)
     ; p34=errorplot(x_grid,dusk_grid(*,6)*(0.75d/Max(dusk_grid(*,6))),err_dusk_grid(*,6)*(0.75d/Max(dusk_grid(*,6))),NAME='Normalized 6716Å',ERRORBAR_COLOR="red",color="red",title='Dusk 5 point Avg 0.04 Resolution',xtitle='Dusk $\rho_c$ ($R_J$)',ytitle='',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,xrange=[5.45,5.85],/overplot,yrange=[0.7,0.83])

     ; legday = LEGEND(TARGET=[p3,p33], POSITION=[5.8,0.75], $
   ;   /DATA, /AUTO_TEXT_COLOR)
   
 



   ;te_fixed_old = interpol(nerneymv1_te,0.01d*dindgen(601) + 4d,x_grid)

   
   

  ; p33=errorplot(x_grid,dusk_grid(*,7),err_dusk_grid(*,7),NAME='6731Å',ERRORBAR_COLOR="brown",color="brown",title='Dusk 5 point Avg 0.04 Resolution',xtitle='Dusk $\rho_c$ ($R_J$)',ytitle='Rayleighs 6731Å',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,xrange=[5.8,6.5],yrange=[0,400])
  ; p2222 = plot(x_grid,dusk_grid(30,7)*(x_grid/x_grid(30))^(-(Ff + Ll)),/overplot)
  ;print, x_grid(30)
   
   ;STOP

  ;legday.save,'5point_Moving_avg_0.04res_dawn_and_dusk_ratios_proper_uncs_4.6_to_7.0.png',resolution=300



  ;stop
  ;;;;;;;;;;;;;;;;;;;;;

  ; indx=where((x_grid ge 4.6d) and (x_grid le 7.0d) );indx=where((x_grid ge 5.) and (x_grid le 6.5) )
  ; x_grid=x_grid(indx) ; values between 5-7 RJ

  ;dawn_grid = dawn_grid(indx,*)

  ;dusk_grid = dusk_grid(indx,*)

  ;err_dawn_grid = err_dawn_grid(indx,*)

  ;err_dusk_grid = err_dusk_grid(indx,*)



  ;err_dawn_grid(*,3) =  5d*err_dawn_grid(*,3)
  ;err_dusk_grid(*,3) =  5d*err_dusk_grid(*,3)
  ;err_dawn_grid(*,4) =  20d*err_dawn_grid(*,4)
  ;err_dusk_grid(*,4) =  20d*err_dusk_grid(*,4)




  p1=errorplot(x_grid,dusk_grid(*,2),ERR_dusk_grid(*,2),NAME='OII (O+) 3729 Angstroms',ERRORBAR_COLOR="blue",color="blue",title='Nominal Coadded APO Dusk Profile From Carl',xtitle='Dusk Distance (RJ)',ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400])

  p2=errorplot(x_grid,dusk_grid(*,6),ERR_dusk_grid(*,6),NAME='SII (S+) 6716 Angstroms',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p3=errorplot(x_grid,dusk_grid(*,7),ERR_dusk_grid(*,7),NAME='SII (S+) 6731 Angstroms',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p4=errorplot(x_grid,dusk_grid(*,0),ERR_dusk_grid(*,0),NAME='SIII (S++) 3722 Angstroms',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p5=errorplot(x_grid,dusk_grid(*,1),ERR_dusk_grid(*,1),NAME='OII (O+) 3726 Angstroms',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p6=errorplot(x_grid,dusk_grid(*,5),ERR_dusk_grid(*,5),NAME='SIII (S++) 6312 Angstroms',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p7=errorplot(x_grid,dusk_grid(*,3),ERR_dusk_grid(*,3),NAME='SII (S+) 4069 Angstroms',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p8=errorplot(x_grid,dusk_grid(*,4),ERR_dusk_grid(*,4),NAME='SII (S+) 4076 Angstroms',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


  leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.1,400.], $
    /DATA, /AUTO_TEXT_COLOR)

  ;stop

  p1=errorplot(x_grid,dawn_grid(*,2),ERR_dawn_grid(*,2),NAME='OII (O+) 3729 Angstroms',ERRORBAR_COLOR="blue",color="blue",title='Nominal Coadded APO Dawn Profile From Carl',xtitle='Dawn Distance (RJ)',ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400])

  p2=errorplot(x_grid,dawn_grid(*,6),ERR_dawn_grid(*,6),NAME='SII (S+) 6716 Angstroms',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p3=errorplot(x_grid,dawn_grid(*,7),ERR_dawn_grid(*,7),NAME='SII (S+) 6731 Angstroms',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p4=errorplot(x_grid,dawn_grid(*,0),ERR_dawn_grid(*,0),NAME='SIII (S++) 3722 Angstroms',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p5=errorplot(x_grid,dawn_grid(*,1),ERR_dawn_grid(*,1),NAME='OII (O+) 3726 Angstroms',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p6=errorplot(x_grid,dawn_grid(*,5),ERR_dawn_grid(*,5),NAME='SIII (S++) 6312 Angstroms',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

  p7=errorplot(x_grid,dawn_grid(*,3),ERR_dawn_grid(*,3),NAME='SII (S+) 4069 Angstroms',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
  p8=errorplot(x_grid,dawn_grid(*,4),ERR_dawn_grid(*,4),NAME='SII (S+) 4076 Angstroms',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


  leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.1,400.], $
    /DATA, /AUTO_TEXT_COLOR)



  ;stop
  ;;;



  idx_max_sp_6731_dusk = where(dusk_grid(*,7) eq max(dusk_grid(*,7)))
  idx_max_sp_6731_dawn = where(dawn_grid(*,7) eq max(dawn_grid(*,7)))

  idx_max_sp_6716_dusk = where(dusk_grid(*,6) eq max(dusk_grid(*,6)))
  idx_max_sp_6716_dawn = where(dawn_grid(*,6) eq max(dawn_grid(*,6)))

  idx_max_s2p_6312_dusk = where(dusk_grid(*,5) eq max(dusk_grid(*,5)))
  idx_max_s2p_6312_dawn = where(dawn_grid(*,5) eq max(dawn_grid(*,5)))

  idx_max_sp_4076_dusk = where(dusk_grid(*,4) eq max(dusk_grid(*,4)))
  idx_max_sp_4076_dawn = where(dawn_grid(*,4) eq max(dawn_grid(*,4)))

  idx_max_sp_4069_dusk = where(dusk_grid(*,3) eq max(dusk_grid(*,3)))
  idx_max_sp_4069_dawn = where(dawn_grid(*,3) eq max(dawn_grid(*,3)))

  idx_max_op_3729_dusk = where(dusk_grid(*,2) eq max(dusk_grid(*,2)))
  idx_max_op_3729_dawn = where(dawn_grid(*,2) eq max(dawn_grid(*,2)))

  idx_max_op_3726_dusk = where(dusk_grid(*,1) eq max(dusk_grid(*,1)))
  idx_max_op_3726_dawn = where(dawn_grid(*,1) eq max(dawn_grid(*,1)))

  idx_max_s2p_3722_dusk = where(dusk_grid(*,0) eq max(dusk_grid(*,0)))
  idx_max_s2p_3722_dawn = where(dawn_grid(*,0) eq max(dawn_grid(*,0)))

  ;print,'$S^+$ 6731 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_sp_6731_dusk),x_grid(idx_max_sp_6731_dawn),x_grid(idx_max_sp_6731_dawn) - x_grid(idx_max_sp_6731_dusk)
  ;print,'$S^+$ 6716 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_sp_6716_dusk),x_grid(idx_max_sp_6716_dawn), x_grid(idx_max_sp_6716_dawn) - x_grid(idx_max_sp_6716_dusk)
  ;print,'$S^+$ 4069 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_sp_4069_dusk),x_grid(idx_max_sp_4069_dawn), x_grid(idx_max_sp_4069_dawn) - x_grid(idx_max_sp_4069_dusk)

  ;print,'$S^+$ 4076 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_sp_4076_dusk),x_grid(idx_max_sp_4076_dawn), x_grid(idx_max_sp_4076_dawn) - x_grid(idx_max_sp_4076_dusk)
  ; print,'$S^{++}$ 6312 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_s2p_6312_dusk),x_grid(idx_max_s2p_6312_dawn), x_grid(idx_max_s2p_6312_dawn) - x_grid(idx_max_s2p_6312_dusk)
  ; print,'$S^{++}$ 3722 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_s2p_3722_dusk),x_grid(idx_max_s2p_3722_dawn), x_grid(idx_max_s2p_3722_dawn) - x_grid(idx_max_s2p_3722_dusk)
  ; print,'$O^{++}$ 3729 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_op_3729_dusk),x_grid(idx_max_op_3729_dawn), x_grid(idx_max_op_3729_dawn) - x_grid(idx_max_op_3729_dusk)
  ; print,'$O^{++}$ 3726 Å peak $\rho_c$ dusk and dawn = ',x_grid(idx_max_op_3726_dusk),x_grid(idx_max_op_3726_dawn), x_grid(idx_max_op_3726_dawn) - x_grid(idx_max_op_3726_dusk)



  ;print,'here'
  ;stop








  ;;;




  ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
  ;x = [3726.,3729.,4069.,6312.,6716.,6731.];

  ;dawn_fit
  norm_vec = [0d,1d,0d] ;looking in yhat direction

  y0 = -7.01d ; just outside outer edge of 6.5d now, 7.0d;7.5d
  z0 = 0d


  x_min=-10d;
  x_max=10d;
  x_step_size= dx_want;0.02d ;  RJ
  n_x=Round((x_max - x_min)/x_step_size) + 1
  xgrid=x_step_size*Dindgen(n_x) + x_min

  y_min=-10d;
  y_max=10d;
  y_step_size=x_step_size ;  RJ
  n_y=Round((y_max - y_min)/y_step_size) + 1
  ygrid=y_step_size*Dindgen(n_y) + y_min

  z_min=-4d;
  z_max=4d;
  z_step_size=x_step_size ;  RJ
  n_z=Round((z_max - z_min)/z_step_size) + 1
  zgrid=z_step_size*Dindgen(n_z) + z_min

  n_x_grid=n_elements(x_grid)


  Te_in=replicate(0d,n_x_grid)


  n_in=replicate({densities, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_x_grid )



  openr,1,'yptsi_8_APO_LINES_optical_vary_nec_tec_lookuptable_66x58x8_CHIANTI_10.1.dat';'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI_10.1.dat' ;  ;'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in=dblarr(66,58,8) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1

  ;nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  ;Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

  nel = [0.000001d,0.00001d,0.0001d,0.001d,0.1d,dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]


  Tec = [0.01,0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]



  ;nelhalf = [replicate(0.5d,9), replicate(5d,9), replicate(50d,4), replicate(125d,39)]

  ;Techalf = [replicate(0.05d,9), replicate(0.25d,17),replicate(0.5d,10),replicate(2.5d,8), replicate(10d,7), replicate(40d,6)]



  n_ne = long64(n_elements(nel))

  n_te = long64(n_elements(tec))

  nel2dto1d =dblarr(n_ne*n_te)
  tec2dto1d =dblarr(n_ne*n_te)
  yptsi_in_2dto1d = dblarr(n_ne*n_te,8)

  ;2nd ->1d mapping for griddata interp
  k=long64(-1)
  for i=0, 66 -1 do begin
    for j=0, 58 -1 do begin
      k = k + long64(1)
      nel2dto1d[k] = nel[i]
      tec2dto1d[k] = tec[j]

      yptsi_in_2dto1d[k,*] = yptsi_in[i,j,*]

      ;print,yptsi_in_2dto1d[k,*]
      ;stop
    endfor
  endfor

  ;write_csv,'nec2dto1d61x57x8.csv',nel2dto1d
  ;write_csv,'tec2dto1d61x57x8.csv',tec2dto1d
  ;write_csv,'ypts2dto1d61x57x8.csv',yptsi_in_2dto1d


  ;stop

  openr,1,'Te_core_general_bagenal1994.csv'
  bag94_te_for_dusk=dblarr(2,133)
  readf,1,bag94_te_for_dusk
  close,1

  print,bag94_te_for_dusk(0,*)
  print,bag94_te_for_dusk(1,*)

  ; p0000  = plot(bag94_te_for_dusk(0,*),bag94_te_for_dusk(1,*),xtitle='Dusk Distance ($R_J$)',ytitle='$T_{ec}$ (eV)', title='Bagenal (1994) Model')
  ;stop
  ;p0000.save,'bag94_Te_used_for_dusk_apo.png',resolution=300

  te_for_plot = dblarr(133)

  for iii=0, 132 do begin

    x000 = bag94_te_for_dusk(0,iii)

      if ( x000 gt 5.7d) then te_for_plot(iii) = interpol(bag94_te_for_dusk(1,*),bag94_te_for_dusk(0,*),x000) else te_for_plot(iii) = 3.79556615d*(x000/5.7d)^(10.1797d) ; for 5.0 to 5.7

      if ( x000 le 5d) then te_for_plot(iii) = 1d

    ;3.79556615d*(x000/5.96d)^(10.6758d) ;for 5.26 to 5.96

    ;if ( x000 gt 5.9d) then te_for_plot(iii) = interpol(bag94_te_for_dusk(1,*) ,bag94_te_for_dusk(0,*) + 0.2d,x000) else te_for_plot(iii) = 3.79556615d*(x000/5.9d)^(10.5614) ; for 5.16 to 5.86

    ;if ( x000 le 5.2d) then te_for_plot(iii) = 1d


  endfor

  x_te_for_plot = reform(bag94_te_for_dusk(0,*))

  x_te_for_plot = [4.5,4.6,4.7,4.8,4.9,x_te_for_plot]
  te_for_plot = [1d,1d,1d,1d,1d,te_for_plot]


  x_grid = reverse(x_grid)
  te_fixed = interpol(te_for_plot,x_te_for_plot,x_grid)
  
  
  openr,1,'Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
  nerneymv1_te=dblarr(601)
  readf,1,nerneymv1_te
  close,1
  
  openr,1,'rT_e_nerney2017.csv'
  rT_e_nerney2017=fltarr(14)
  readf,1,rT_e_nerney2017
  close,1
  
  openr,1,'T_e_nerney2017.csv'
  T_e_nerney2017=fltarr(14)
  readf,1,T_e_nerney2017
  close,1
  
  openr,1,'rsteffl2004b_Te.csv'
  rT_e_steffl2004b=fltarr(12)
  readf,1,rT_e_steffl2004b
  close,1

  openr,1,'steffl2004b_Te.csv'
  T_e_steffl2004b=fltarr(12)
  readf,1,T_e_steffl2004b
  close,1




  
  openr,1,'cm3_proper_runningavg5points_dl=0.04_p_output_dusk_only_high_res_lines_limitTe_le8_ftol-20_into_4.6_out_to7.csv'
  p_dusky=fltarr(60,5)
  readf,1,p_dusky
  close,1
  
  
  
  te_nerneymv1_interped = interpol(nerneymv1_te,0.01d*dindgen(601) + 4d,x_grid)
  
  te_steffl2004b_interped = interpol( T_e_steffl2004b,rT_e_steffl2004b,x_grid)
  te_nerney2017_interped = interpol( T_e_nerney2017,rT_e_nerney2017,x_grid)

  te_fixed =  te_steffl2004b_interped;te_nerney2017_interped;
   te_fixed(21:59) = te_nerneymv1_interped(21:59)
   
   
   for iii=0, n_elements(x_grid) -1 do begin
     if ( x_grid(iii) le 5.9d) then te_fixed(iii) = 3.79556615d*(x_grid(iii)/5.7d)^(10.1797d) ; for 5.24 to 5.94
     if ( x_grid(iii) le 5.0d) then te_fixed(iii) = 1d
   endfor
   ;te_fixed = 1.2d*te_fixed
  
  ;te_fixed = p_dusky(*,0)
 
 ;p1=plot(x_grid,te_fixed)
 ; p1=plot(x_grid,3.79556615d*(x_grid/5.7d)^(10.1797d),/overplot,color='blue' )
 ; stop

  ;  p00000  = plot(x_te_for_plot,te_for_plot,xtitle='Dawn Distance ($R_J$)',ytitle='$T_{ec}$ (eV)', title='Modified for APO Dawn (0.16 $R_J$ Offset) use from Bagenal (1994) Model',xrange=[4.5,7.0],yrange=[0,8])
  ;stop
  ;p00000.save,'Modified_for_APO_Dawn_strict_0.16_offset_use_from_bag94_Tec.png',resolution=300

  ;stop
  ;
  
  
  
  x_grid_coadded = x_grid
  ;xgrid_strings = ['x_grid_night1.csv','x_grid_night2.csv','x_grid_night3.csv','x_grid_night4.csv','x_grid_dawn_night5.csv','x_grid_night6.csv']
  ;xgrid_strings = ['x_grid_night1.csv','xgrid_night2_fixed.csv','x_grid_night3.csv','x_grid_night4.csv','x_grid_dawn_night5.csv','x_grid_night6.csv']
  xgrid_strings = ['xgrid_night1_fixed.csv','xgrid_night2_fixed.csv','xgrid_night3_fixed.csv','xgrid_night4_fixed.csv','xgrid_night5_dawn_fixed.csv','xgrid_night6_fixed.csv']
  nightname_strings = ['night1','night2','night3','night4','night5','night6']
  spec_strings = ['Dawn_night1_grid_to_plot_71x8.csv','Dawn_night2_grid_to_plot_71x8.csv','Dawn_night3_grid_to_plot_68x8.csv','Dawn_night4_grid_to_plot_67x8.csv','Dawn_night5_grid_to_plot_55x8.csv','Dawn_night6_grid_to_plot_64x8.csv']

  err_spec_strings =  ['err_Dawn_night1_grid_to_plot_71x8.csv','err_Dawn_night2_grid_to_plot_71x8.csv','err_Dawn_night3_grid_to_plot_68x8.csv','err_Dawn_night4_grid_to_plot_67x8.csv','err_Dawn_night5_grid_to_plot_55x8.csv','err_Dawn_night6_grid_to_plot_64x8.csv']

  infile_strings = ['night1_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv','night2_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv',$
    'night3_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv','night4_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv', $
    'night5_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv','night6_varypls_varytec_locs_varyrhoRnsp_and_rhoRns2p_and_nop_diff_ev_outs_and_errors_dawn_funcform_tec3partpl_python_60x16.csv']
  nxxx = [71,71,68,67,55,64]

  for iiii = 4, 4 do begin;n_elements(nxxx) - 1 do begin

;iiii = 1;0

    current_night = nightname_strings[iiii]
    current_n_xgrid = nxxx[iiii]
    x_grid_string = xgrid_strings[iiii]
    spec_string = spec_strings[iiii]
    err_spec_string = err_spec_strings[iiii]

    infile_string = infile_strings[iiii]

    openr,1,x_grid_string ;'x_grid_dusk_night5.csv'
    x_grid = dblarr(current_n_xgrid)
    readf,1,x_grid
    close,1
    
    x_step_size = abs(x_grid(0) - x_grid(1))

    openr,1,spec_string ;'x_grid_dusk_night5.csv'
    dawn_grid = dblarr(current_n_xgrid,8)
    readf,1,dawn_grid
    close,1

    openr,1,err_spec_string ;'x_grid_dusk_night5.csv'
    err_dawn_grid = dblarr(current_n_xgrid,8)
    readf,1,err_dawn_grid
    close,1

  ;dawn_grid = dusk_grid
  ;err_dawn_grid = err_dusk_grid
   XTITTLEE = 'Dawn $\rho_c$ ($R_J$)'
    idx_want= [0,1,2,3,5,6,7];[3,6,7]
    dawn_grid_temp = dawn_grid(*,idx_want);so currently 7 lines ;
    err_dawn_grid_temp = err_dawn_grid(*,idx_want);;so currently 7 lines (*,idx_want)


    ;y needs to be setup to for [model] =  reform(dblarr(num_elements,8),num_elements*8)
    ndawn_grid = n_elements(dawn_grid_temp)
    y=reform(dawn_grid_temp,ndawn_grid);[dawn_grid(n_elements(x_grid)-1-i,3),dawn_grid(n_elements(x_grid)-1-i,6),dawn_grid(n_elements(x_grid)-1-i,7)]
    err=reform(err_dawn_grid_temp,ndawn_grid);[50d*err_dawn_grid(n_elements(x_grid)-1-i,3),err_dawn_grid(n_elements(x_grid)-1-i,6),err_dawn_grid(n_elements(x_grid)-1-i,7)]
    x = replicate(1d,ndawn_grid)

openr,1,'dawn_apo_final_totchmixratio.csv'
totchmix_for_apo_new_final_in = dblarr(60) ;
readf,1,totchmix_for_apo_new_final_in
close,1

totchmix_for_apo = totchmix_for_apo_new_final_in

totchmix_for_apo = interpol(totchmix_for_apo_new_final_in,x_grid_coadded,x_grid)
totchmix_for_apo[-1] = 1d

openr,1,infile_string
outin = dblarr(16,current_n_xgrid) ;
readf,1,outin
close,1
tec = reform(outin[0,*])
nec_funcform_test = reform(outin[2,*])
nsp_funcform_test = reform(outin[4,*])
ns2p_funcform_test = reform(outin[6,*])
nop_funcform_test = reform(outin[8,*])






x_grid_higher_res = x_grid
   
model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(tec, nec_funcform_test,nsp_funcform_test,ns2p_funcform_test,nop_funcform_test, max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
  

  yptsi_output1 = model
    
   ; p11=errorplot(x_grid,dawn_grid(*,2),ERR_dawn_grid(*,2),NAME='OII (O+) 3729 Å',ERRORBAR_COLOR="blue",color="blue",xtitle=  XTITTLEE,ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,xrange=[4.5,7.0],yrange=[0,400],ASPECT_RATIO=0.0093975)
   ; p111 = plot(x_grid_higher_res,yptsi_output1(*,2),color="blue",/overplot)
   ; p22=errorplot(x_grid,dawn_grid(*,6),ERR_dawn_grid(*,6),NAME='SII (S+) 6716 Å',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p222 = plot(x_grid_higher_res,yptsi_output1(*,6),color="red",/overplot)
   ; p33=errorplot(x_grid,dawn_grid(*,7),ERR_dawn_grid(*,7),NAME='SII (S+) 6731 Å',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p333 = plot(x_grid_higher_res,yptsi_output1(*,7),color="brown",/overplot)
   ; p44=errorplot(x_grid,dawn_grid(*,0),ERR_dawn_grid(*,0),NAME='SIII (S++) 3722 Å',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p444 = plot(x_grid_higher_res,yptsi_output1(*,0),color="pink",/overplot)
   ; p55=errorplot(x_grid,dawn_grid(*,1),ERR_dawn_grid(*,1),NAME='OII (O+) 3726 Å',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p555 = plot(x_grid_higher_res,yptsi_output1(*,1),color="purple",/overplot)
   ; p66=errorplot(x_grid,dawn_grid(*,5),ERR_dawn_grid(*,5),NAME='SIII (S++) 6312 Å',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p666 = plot(x_grid_higher_res,yptsi_output1(*,5),color="orange",/overplot)
   ; p77=errorplot(x_grid,dawn_grid(*,3),ERR_dawn_grid(*,3),NAME='SII (S+) 4069 Å',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p777 = plot(x_grid_higher_res,yptsi_output1(*,3),color="light Green",/overplot)
   ; p88=errorplot(x_grid,dawn_grid(*,4),ERR_dawn_grid(*,4),NAME='SII (S+) 4076 Å',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
   ; p888 = plot(x_grid_higher_res,yptsi_output1(*,4),color="dark green",/overplot)

   ; leg22 = LEGEND(TARGET=[p33,p22,p66,p88,p77,p11,p55,p44], POSITION=[7.2,400.], $
   ;   /DATA, /AUTO_TEXT_COLOR)

    ;leg22.save, 'dusk_coadded_plug_in_smoothed_nec_nsp1_test_drho=0.01_onion_sum_method.png', RESOLUTION=300
   
   tec_fixed_found = tec 
   nec_fixed_found =  nec_funcform_test 
    nsp_fixed_found = nsp_funcform_test 
   ns2p_fixed_found =  ns2p_funcform_test 
   nop_fixed_found =  nop_funcform_test 
   
   nspmix =  nsp_fixed_found/nec_fixed_found
   ns2pmix =  ns2p_fixed_found/nec_fixed_found
   nopmix =  nop_fixed_found/nec_fixed_found
   
     nec_fixed_found_og = nec_fixed_found 
     
     tec_fixed_found_og = tec_fixed_found
     
     nsp_fixed_found_og = nsp_funcform_test
     ns2p_fixed_found_og =  ns2p_funcform_test
     nop_fixed_found_og =  nop_funcform_test
     
     nspmix_og =  nsp_fixed_found/nec_fixed_found
     ns2pmix_og =  ns2p_fixed_found/nec_fixed_found
     nopmix_og =  nop_fixed_found/nec_fixed_found
     
     left = 0.2
     right = 0.05
     MARGINz0 =[left,0.1,right,0.1]
     MARGINz1 =[left,0.1,right,0.0]
     MARGINz2 =[left,0.25,right,0.0]
     
    ; p1=plot(x_grid,tec_fixed_found,ytitle='$T_{ec}$ (eV)',layout=[2,4,1],margin =   MARGINz0,XTICKFORMAT="(A1)",yrange=[0,10]);,/ylog)
    ; p2=plot(x_grid,nec_fixed_found,ytitle='$n_{ec}$ ($cm^{-3}$)',layout=[2,4,2],/current,margin =   MARGINz0,XTICKFORMAT="(A1)",yrange=[500,4000]);,/ylog)
    ; p3=plot(x_grid,nsp_fixed_found,ytitle='$n_{S^{+}}$ ($cm^{-3}$)',layout=[2,4,3],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[0,800]);,/ylog)
    ; p4=plot(x_grid,nsp_fixed_found/nec_fixed_found,ytitle='$n_{S^{+}}/n_{ec}$ ',layout=[2,4,4],/current,margin =   MARGINz1,yrange=[0,0.8],XTICKFORMAT="(A1)");,/ylog)
    ; p5=plot(x_grid,ns2p_fixed_found,ytitle='$n_{S^{++}}$ ($cm^{-3}$)',layout=[2,4,5],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[0,800]);,/ylog)
    ; p6=plot(x_grid,ns2p_fixed_found/nec_fixed_found,ytitle='$n_{S^{++}}/n_{ec}$ ',layout=[2,4,6],/current,margin =   MARGINz1,yrange=[0,0.4],XTICKFORMAT="(A1)");,/ylog)
    ; p7=plot(x_grid,nop_fixed_found,ytitle='$n_{O^{+}}$ ($cm^{-3}$)',layout=[2,4,7],/current,margin =   MARGINz2,yrange=[0,800],xtitle=   XTITTLEE);,/ylog)
    ; p8=plot(x_grid,nop_fixed_found/nec_fixed_found,xtitle=   XTITTLEE,ytitle='$n_{O^{+}}/n_{ec}$ ',layout=[2,4,8],/current,margin =   MARGINz2,yrange=[0,0.6]);,/ylog)

   xtest = double([6.982  , 6.94828, 6.91456, 6.88084, 6.84712, 6.8134 , 6.77968, $
     6.74596, 6.71224, 6.67852, 6.6448 , 6.61108, 6.57736, 6.54364, $
     6.50992, 6.4762 , 6.44248, 6.40876, 6.37504, 6.34132, 6.3076 ,$
     6.27388, 6.24016, 6.20644, 6.17272, 6.139  , 6.10528, 6.07156,$
     6.03784, 6.00412, 5.9704 , 5.93668, 5.90296, 5.86924, 5.83552,$
     5.8018 , 5.76808, 5.73436, 5.70064, 5.66692, 5.6332 , 5.59948,$
     5.56576, 5.53204, 5.49832, 5.4646 , 5.43088, 5.39716, 5.36344,$
     5.32972, 5.296  , 5.26228, 5.22856, 5.19484, 5.16112, 5.1274 ,$
     5.09368, 5.05996, 5.02624, 4.99252, 4.9588 , 4.92508, 4.89136,$
     4.85764, 4.82392, 4.7902 , 4.75648, 4.72276, 4.68904, 4.65532,$
     4.6216 ])
     
     nec_testing = double([1309.52839713, 1321.60903471, 1333.81142346, 1346.15093037,$
   1358.66316191, 1371.38944938, 1384.36993771, 1397.6115489 ,$
   1411.12032573, 1424.90278621, 1438.96569753, 1453.31603904,$
   1467.96099968, 1482.9089847 , 1497.31746246, 1511.15397385,$
   1524.42043433, 1536.29759826, 1546.11590063, 1556.2467663 ,$
   1566.6947353 , 1577.46459408, 1588.56138002, 1599.9903863 ,$
   1611.75716733, 1623.86754447, 1646.6877167 , 1679.51778618,$
   1703.05533717, 1715.92466145, 1717.00834813, 1705.60315376,$
   1681.63157042, 1645.81770899, 1593.75395316, 1520.54834668,$
   1434.7206971 , 1345.20601609, 1262.88439453, 1196.02485931,$
   1149.32418856, 1124.36283325, 1120.59899068, 1136.27051398,$
   1168.89210846, 1215.36589936, 1271.90008048, 1330.79151055,$
   1386.65877618, 1435.29987793, 1472.98845331, 1494.94797805,$
   1497.17577199, 1477.09875831, 1431.96699695, 1356.61860349,$
   1255.28118348, 1133.6734106 ,  997.12011015,  853.51679665,$
   710.75022712,  575.62717625,  453.31281235,  347.07833735,$
   258.34382672,  186.93178819,  131.48469151,   89.90571813,$
   59.76544712,   38.62913901,   24.27759267])
   
   nsp_testing = double([ 74.89787847,  76.06617409,  77.2585106 ,  78.47550359,$
 79.71778754,  80.98601649,  82.28086476,  83.60302767,$
 84.95322229,  86.33218827,  87.74068865,  89.17951072,$
 90.64946691,  92.15139576,  93.68616288,  95.25466195,$
 96.85781582,  98.49657761, 100.17193186, 101.88489572,$
 103.63652026, 105.42789177, 107.26013311, 109.13440519,$
 111.05190845, 113.01388444, 123.85669564, 148.6526711 ,$
 173.68892184, 197.62450264, 219.06251723, 236.72583291,$
 249.64479884, 257.32792253, 254.95396555, 237.05257533,$
 209.75403403, 181.25379695, 159.30831076, 149.45607828,$
 154.38451048, 174.31720461, 207.89632294, 253.02300723,$
 307.36352177, 368.5044409 , 433.89923325, 500.77184429,$
 566.08791898, 626.63468007, 679.19989737, 720.81419162,$
 749.01195944, 762.06697656, 757.78691912, 732.33229074,$
 687.54789803, 627.09224015, 555.64010874, 478.28795222,$
 399.96223632, 324.9243265 , 256.43609785, 196.61189304,$
 146.4448497 , 105.96735966,  74.49100213,  50.87087565,$
 33.74957455,  21.75209466,  13.61970056])
   
   ns2p_testing = double([2.65637155e+02, 2.69190173e+02, 2.72808346e+02, 2.76493198e+02,$
 2.80246291e+02, 2.84069237e+02, 2.87963690e+02, 2.91931353e+02,$
 2.95973979e+02, 3.00093372e+02, 3.04291388e+02, 3.08569939e+02,$
 3.12930991e+02, 3.17376571e+02, 3.21908766e+02, 3.26529725e+02,$
 3.31241662e+02, 3.36046859e+02, 3.40947667e+02, 3.45946510e+02,$
 3.51045886e+02, 3.56248369e+02, 3.61556617e+02, 3.66973368e+02,$
 3.72501447e+02, 3.78143771e+02, 3.82294264e+02, 3.82032014e+02,$
 3.78684809e+02, 3.72312840e+02, 3.63069819e+02, 3.51175564e+02,$
 3.36907309e+02, 3.20588974e+02, 3.02578966e+02, 2.83257157e+02,$
 2.63011692e+02, 2.42226259e+02, 2.21268400e+02, 2.00479316e+02,$
 1.80165550e+02, 1.60592751e+02, 1.41981621e+02, 1.24506013e+02,$
 1.08293035e+02, 9.34249189e+01, 7.99423552e+01, 6.78489534e+01,$
 5.71164717e+01, 4.76904691e+01, 3.94960644e+01, 3.24435265e+01,$
 2.64334770e+01, 2.13615388e+01, 1.71223268e+01, 1.36127238e+01,$
 1.07344374e+01, 8.39586696e+00, 6.51334096e+00, 5.01180405e+00,$
 3.82504332e+00, 2.89554712e+00, 2.17408701e+00, 1.61910573e+00,$
 1.19598436e+00, 8.76249533e-01, 6.36769278e-01, 4.58974343e-01,$
 3.28130591e-01, 2.32678919e-01, 1.63651262e-01])
   
   nop_testing = double([243.93171785, 247.64155468, 251.42627406, 255.28776381,$
 259.22796889, 263.24889348, 267.35260306, 271.54122661,$
 275.81695888, 280.18206278, 284.63887185, 289.18979284,$
 293.83730839, 298.58397984, 303.43245014, 308.38544689,$
 313.4457855 , 318.61637254, 323.90020914, 329.30039461,$
 334.8201302 , 340.46272302, 346.23159014, 352.13026284,$
 358.16239113, 364.33174836, 372.96143941, 384.48403133,$
 395.11631921, 404.82218066, 413.60393403, 421.51398314,$
 428.66771223, 435.25663968, 441.5601607 , 447.95346171,$
 454.19437326, 458.54747439, 461.61310968, 464.32279614,$
 467.68468401, 472.68813079, 480.1929624 , 490.81367423,$
 504.81200056, 522.0130833 , 541.76017284, 562.91992534,$
 583.94489372, 602.99227318, 618.08943443, 627.32873582,$
 629.06813846, 622.11160491, 604.86297192, 574.4977692 ,$
 532.19962197, 480.41567671, 422.30031007, 361.30173227,$
 300.74890033, 243.50527784, 191.73523658, 146.80151867,$
 109.28546387,  79.10161329,  55.66803478,  38.09340303,$
 25.34918888,  16.4066983 ,  10.33058959])
 
 tec_testing = double([4.05989185, 4.02476108, 3.98976562, 3.95490559, 3.92018114,$
 3.8855924 , 3.8511395 , 3.81682259, 3.78264181, 3.74859729,$
 3.71468918, 3.68091762, 3.64728274, 3.6137847 , 3.58042363,$
 3.54719969, 3.51411302, 3.48116376, 3.44835206, 3.41567807,$
 3.38314195, 3.35074383, 3.31848388, 3.28636225, 3.25437909,$
 3.22253455, 3.19082879, 3.15926196, 3.12783423, 3.09654576,$
 2.88873892, 2.66083243, 2.44975953, 2.25436254, 2.0735574 ,$
 1.90632944, 1.75172927, 1.60886891, 1.47691817, 1.35510117,$
 1.24269312, 1.16365757, 1.10153828, 1.04238752, 0.98608024,$
 0.93249593, 0.88151852, 0.83303619, 0.78694128, 0.74313012,$
 0.70150296, 0.6619638 , 0.62442029, 0.5887836 , 0.55496834,$
 0.52289239, 0.49247688, 0.46364597, 0.43632687, 0.41044966,$
 0.38594719, 0.36275505, 0.34081143, 0.32005703, 0.30043499,$
 0.2818908 , 0.26437223, 0.24782923, 0.23221385, 0.2174802 ,$
 0.20358434])
     print,max(abs(x_grid - xtest))
       print,max(abs(tec_fixed_found - tec_testing))
      print,max(abs(nec_fixed_found - nec_testing))
      print,max(abs(nsp_fixed_found - nsp_testing))
      print,max(abs(ns2p_fixed_found - ns2p_testing))
      print,max(abs(nop_fixed_found - nop_testing))
  ; stop
    
    ;what used before to start
    ;lambda1_arr = [1d-7];[1d-5,0.2d-6,0.4d-6,0.6d-6,0.8d-6,1d-6,2d-7,4d-7,6d-7,8d-7,1d-7] ;[1d-6,1d-8,1d-10,1d-12,1d-7,1d-9,1d-11]

   ; lambda2_arr = [1d-2];[2d-2,4d-2,6d-2,8d-2,1d-2,2d-3,4d-3,6d-3,8d-3,1d-3];[1d-4] ;[1d-1,5d-2,1d-2,5d-3,1d-3]


    ;lambda1_arr = [5d-5,1d-5,5d-6,1d-6,5d-7,1d-7,5d-8];[1000d,100d,10d,1d,0.1,0.01d,0.001d,0.0001d]; [1d-6,1d-7,1d-8,1d-9,1d-12,1d-18,1d-24];[1d-5,0.2d-6,0.4d-6,0.6d-6,0.8d-6,1d-6,2d-7,4d-7,6d-7,8d-7,1d-7] ;[1d-6,1d-8,1d-10,1d-12,1d-7,1d-9,1d-11]

    ;lambda2_arr = [1d-1,5d-2,1d-2,5d-3,1d-3,5d-4,1d-4];[1000d,100d,10d,1d,0.1,0.01d,0.001d,0.0001d];[1d-3,1d-4,1d-5,1d-6,1d-7,1d-8,1d-9,1d-10,1d-15,1d-18,1d-24];[2d-2,4d-2,6d-2,8d-2,1d-2,2d-3,4d-3,6d-3,8d-3,1d-3];[1d-4] ;[1d-1,5d-2,1d-2,5d-3,1d-3]

    lambda1_arr = [1d-7,1d-8];[1d-8,1d-9,1d-10];[5d-5,1d-5,5d-6,1d-6,5d-7,1d-7,5d-8];[1000d,100d,10d,1d,0.1,0.01d,0.001d,0.0001d]; [1d-6,1d-7,1d-8,1d-9,1d-12,1d-18,1d-24];[1d-5,0.2d-6,0.4d-6,0.6d-6,0.8d-6,1d-6,2d-7,4d-7,6d-7,8d-7,1d-7] ;[1d-6,1d-8,1d-10,1d-12,1d-7,1d-9,1d-11]

    lambda2_arr = [5d-5,1d-5,5d-6,1d-7];[1d-2,5d-3,1d-3,1d-4];[1d-1,5d-2,1d-2,5d-3,1d-3,5d-4,1d-4];[1000d,100d,10d,1d,0.1,0.01d,0.001d,0.0001d];[1d-3,1d-4,1d-5,1d-6,1d-7,1d-8,1d-9,1d-10,1d-15,1d-18,1d-24];[2d-2,4d-2,6d-2,8d-2,1d-2,2d-3,4d-3,6d-3,8d-3,1d-3];[1d-4] ;[1d-1,5d-2,1d-2,5d-3,1d-3]


    y_old = y
    err_old = err
    x_old = x
   
   for iiiii = 0, n_elements(lambda1_arr) - 1 do begin
    for jjjjj = 0, n_elements(lambda2_arr) - 1 do begin

    
 lambdaa1 = lambda1_arr[iiiii]
 
 lambdaa2 =   lambda2_arr[jjjjj];lambdaa1;*lambdafac 
   
   num_of_fit_params = 2
   n_x_grid = n_elements(x_grid)
   nx_grid = n_x_grid
   parinfo = replicate({value:0.d, fixed:0, limited:[0,0], limits:[0.,0]}, num_of_fit_params*n_x_grid);;;!!????????? below fix

  

   parinfo[*].value = [nec_fixed_found_og,tec_fixed_found_og]


   parinfo[0:1*nx_grid - 1].limited[0] = 1
   parinfo[0:1*nx_grid - 1].limits[0]  = 0.00001d;0.1;d
   parinfo[0:1*nx_grid - 1].limited[1] = 1
   parinfo[0:1*nx_grid - 1].limits[1]  = 5000.d;d

   parinfo[1*nx_grid:2*nx_grid - 1].limited[0] = 1
   parinfo[1*nx_grid:2*nx_grid - 1].limits[0]  = 0.1d;0.1;d
   parinfo[1*nx_grid:2*nx_grid - 1].limited[1] = 1
   parinfo[1*nx_grid:2*nx_grid - 1].limits[1]  = 10.d;d
   
   
   y=[y_old,replicate(1d,2*n_elements(x_grid))]
   err=[err_old,replicate(1d,2*n_elements(x_grid))]
   x = [x_old,replicate(1d,2*n_elements(x_grid))]
   fa = {X:x, Y:y, ERR:err}
   
   p=mpfit('apo_deviates_CSR_fixed_mixr',functargs=fa,parinfo=parinfo,perror=perror);,MAXITER=50);,ftol=1e-20,xtol=1e-20,gtol=1e-20);,MAXITER=0);,ftol=1e-40,xtol=1e-40,gtol=1e-40,MAXITER=0)
   

   
   nec_fixed_found = p[0:nx_grid - 1]
   
   tec_fixed_found = p[1*nx_grid:2*nx_grid - 1]
   
   nsp_fixed_found = nspmix * nec_fixed_found
   ns2p_fixed_found = ns2pmix * nec_fixed_found
   nop_fixed_found = nopmix * nec_fixed_found

   y = y_old
   err = err_old
   x = x_old



num_of_fit_params = 4
n_x_grid = n_elements(x_grid)
nx_grid = n_x_grid
parinfo = replicate({value:0.d, fixed:0, limited:[0,0], limits:[0.,0]}, num_of_fit_params*n_x_grid);;;!!????????? below fix


parinfo[*].value = [nsp_fixed_found,ns2p_fixed_found,nop_fixed_found,tec_fixed_found]


parinfo[0:3*nx_grid - 1].limited[0] = 1
parinfo[0:3*nx_grid - 1].limits[0]  = 0.d;0.1;d
parinfo[0:3*nx_grid - 1].limited[1] = 1
parinfo[0:3*nx_grid - 1].limits[1]  = 5000.d;d

parinfo[3*nx_grid:4*nx_grid - 1].limited[0] = 1
parinfo[3*nx_grid:4*nx_grid - 1].limits[0]  = 0.1d;0.1;d
parinfo[3*nx_grid:4*nx_grid - 1].limited[1] = 1
parinfo[3*nx_grid:4*nx_grid - 1].limits[1]  = 10.d;d

fa = {X:x, Y:y, ERR:err}

p=mpfit('apo_deviates_emission_table_allside_simult_foursix_to_seven_7lines_no_func_form_fixed_necbyratio_tec_pl_double_dawn',functargs=fa,parinfo=parinfo,perror=perror,MAXITER=0);,ftol=1e-20,xtol=1e-20,gtol=1e-20);,MAXITER=0);,ftol=1e-40,xtol=1e-40,gtol=1e-40,MAXITER=0)

nsp_fixed_found = p[0:nx_grid - 1]
ns2p_fixed_found = p[nx_grid:2*nx_grid - 1]
nop_fixed_found = p[2*nx_grid:3*nx_grid - 1]
nsp_error = perror[0:nx_grid - 1]
ns2p_error = perror[nx_grid:2*nx_grid - 1]
nop_error = perror[2*nx_grid:3*nx_grid - 1]

tec_fixed_found = p[3*nx_grid:4*nx_grid - 1]
tec_error = perror[3*nx_grid:4*nx_grid - 1]



nec_fixed_found = (nsp_fixed_found  + 2d*ns2p_fixed_found  + nop_fixed_found)/ totchmix_for_apo 


nec_error = sqrt((nsp_error^2d) + 4d*(ns2p_error^2d) + (nop_error^2d)) / totchmix_for_apo 

nspmix_error = (nsp_fixed_found/nec_fixed_found)*Sqrt((nsp_error/nsp_fixed_found)^2d + (nec_error/nec_fixed_found)^2d)

ns2pmix_error = (ns2p_fixed_found/nec_fixed_found)*Sqrt((ns2p_error/ns2p_fixed_found)^2d + (nec_error/nec_fixed_found)^2d)

nopmix_error = (nop_fixed_found/nec_fixed_found)*Sqrt((nop_error/nop_fixed_found)^2d + (nec_error/nec_fixed_found)^2d)




totchmix = (nsp_fixed_found + 2d*ns2p_fixed_found + nop_fixed_found) / nec_fixed_found

totchmix_error =  totchmix*Sqrt( (nsp_error^2d + 4d*ns2p_error^2d + nop_error^2d)/(nsp_fixed_found + 2d*ns2p_fixed_found + nop_fixed_found)^2d + (nec_error/nec_fixed_found)^2d)


model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(tec_fixed_found, nec_fixed_found,nsp_fixed_found,ns2p_fixed_found,nop_fixed_found, max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
yptsi_output1 = model
p11=errorplot(x_grid,dawn_grid(*,2),ERR_dawn_grid(*,2),NAME='OII (O+) 3729 Å',ERRORBAR_COLOR="blue",color="blue",xtitle=  XTITTLEE,ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400])
p111 = plot(x_grid_higher_res,yptsi_output1(*,2),color="blue",/overplot)
p22=errorplot(x_grid,dawn_grid(*,6),ERR_dawn_grid(*,6),NAME='SII (S+) 6716 Å',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p222 = plot(x_grid_higher_res,yptsi_output1(*,6),color="red",/overplot)
p33=errorplot(x_grid,dawn_grid(*,7),ERR_dawn_grid(*,7),NAME='SII (S+) 6731 Å',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p333 = plot(x_grid_higher_res,yptsi_output1(*,7),color="brown",/overplot)
p44=errorplot(x_grid,dawn_grid(*,0),ERR_dawn_grid(*,0),NAME='SIII (S++) 3722 Å',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p444 = plot(x_grid_higher_res,yptsi_output1(*,0),color="pink",/overplot)
p55=errorplot(x_grid,dawn_grid(*,1),ERR_dawn_grid(*,1),NAME='OII (O+) 3726 Å',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p555 = plot(x_grid_higher_res,yptsi_output1(*,1),color="purple",/overplot)
p66=errorplot(x_grid,dawn_grid(*,5),ERR_dawn_grid(*,5),NAME='SIII (S++) 6312 Å',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p666 = plot(x_grid_higher_res,yptsi_output1(*,5),color="orange",/overplot)
p77=errorplot(x_grid,dawn_grid(*,3),ERR_dawn_grid(*,3),NAME='SII (S+) 4069 Å',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p777 = plot(x_grid_higher_res,yptsi_output1(*,3),color="light Green",/overplot)
p88=errorplot(x_grid,dawn_grid(*,4),ERR_dawn_grid(*,4),NAME='SII (S+) 4076 Å',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p888 = plot(x_grid_higher_res,yptsi_output1(*,4),color="dark green",/overplot)

leg22 = LEGEND(TARGET=[p33,p22,p66,p88,p77,p11,p55,p44], POSITION=[7.2,400.], $
  /DATA, /AUTO_TEXT_COLOR)


formatted_lambdaa1 = string(lambdaa1, FORMAT='(E10.2)')
trimmed_lambdaa1 = strtrim(formatted_lambdaa1, 1)
formatted_lambdaa2 = string(lambdaa2, FORMAT='(E10.2)')
trimmed_lambdaa2 = strtrim(formatted_lambdaa2, 1)

leg22.save,current_night + '_lambda1=' + trimmed_lambdaa1  + '_lambda2=' + trimmed_lambdaa2 + '_CSR_dawn_diffev_python_model_vs_obs.png',resolution=300

;doubles_L2_regularization_lambdaeq0.00017_let_it
;0it


left = 0.2
right = 0.05
MARGINz0 =[left,0.1,right,0.1]
MARGINz1 =[left,0.1,right,0.0]
MARGINz2 =[left,0.25,right,0.0]


p1=errorplot(x_grid,tec_fixed_found,tec_error,ytitle='$T_{ec}$ (eV)',layout=[2,4,1],margin =   MARGINz0,XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,10]);,/ylog)
p2=errorplot(x_grid,nec_fixed_found,nec_error,ytitle='$n_{ec}$ ($cm^{-3}$)',layout=[2,4,2],/current,margin =   MARGINz0,XTICKFORMAT="(A1)",yrange=[500,4000],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p3=errorplot(x_grid,nsp_fixed_found,nsp_error,ytitle='$n_{S^{+}}$ ($cm^{-3}$)',layout=[2,4,3],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[0,800],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p4=errorplot(x_grid,nsp_fixed_found/nec_fixed_found,nspmix_error,ytitle='$n_{S^{+}}/n_{ec}$ ',layout=[2,4,4],/current,margin =   MARGINz1,yrange=[0,0.8],XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p5=errorplot(x_grid,ns2p_fixed_found,ns2p_error,ytitle='$n_{S^{++}}$ ($cm^{-3}$)',layout=[2,4,5],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[0,800],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p6=errorplot(x_grid,ns2p_fixed_found/nec_fixed_found,ns2pmix_error,ytitle='$n_{S^{++}}/n_{ec}$ ',layout=[2,4,6],/current,margin =   MARGINz1,yrange=[0,0.4],XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p7=errorplot(x_grid,nop_fixed_found,nop_error,ytitle='$n_{O^{+}}$ ($cm^{-3}$)',layout=[2,4,7],/current,margin =   MARGINz2,yrange=[0,800],xtitle=   XTITTLEE,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p8=errorplot(x_grid,nop_fixed_found/nec_fixed_found,nopmix_error,xtitle=   XTITTLEE,ytitle='$n_{O^{+}}/n_{ec}$ ',layout=[2,4,8],/current,margin =   MARGINz2,yrange=[0,0.6],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)

;print,min(te_fixed),min(nec_fixed_found),min(nsp_fixed_found),min(ns2p_fixed_found),min(nop_fixed_found)


p8.save,current_night + '_lambda1=' + trimmed_lambdaa1  + '_lambda2=' + trimmed_lambdaa2 + '_CSR_dawn_diffev_python_simult_param_errorplots.png',resolution=300

output = [[tec_fixed_found],[tec_error],[nec_fixed_found],[nec_error],[nsp_fixed_found],[nsp_error],[ns2p_fixed_found],[ns2p_error],[nop_fixed_found],[nop_error],[nspmix_error],[ns2pmix_error],[nopmix_error]]

write_csv,current_night + '_lambda1=' + trimmed_lambdaa1  + '_lambda2=' + trimmed_lambdaa2 + '_CSR_dawn_diffev_output_to_plot'+ strtrim(string(current_n_xgrid),1)+'x13.csv',output

write_csv,current_night + '_lambda1=' + trimmed_lambdaa1  + '_lambda2=' + trimmed_lambdaa2 + '_CSR_dawn_diffev_model_to_plot'+ strtrim(string(current_n_xgrid),1)+'x8.csv',model



save,tec_fixed_found,tec_error,nec_fixed_found,nec_error,nsp_fixed_found,nsp_error,ns2p_fixed_found,ns2p_error,nop_fixed_found,nop_error,x_grid,nspmix_error,ns2pmix_error,nopmix_error,totchmix_for_apo,lambdaa1,lambdaa2, FILENAME = current_night + '_lambda1=' + trimmed_lambdaa1  + '_lambda2=' + trimmed_lambdaa2 + '_CSR_dawn_diffev_python_simult_param_save_file_simult_errors_from_mpfit_no_polish.sav',/verbose
;stop
 endfor
    endfor
    
    ;stop
    endfor

 stop


end



