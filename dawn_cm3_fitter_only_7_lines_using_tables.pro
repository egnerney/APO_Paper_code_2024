


Function APO_cm3_model_given_emiss_table,Tec_a,nec_a,Nsp_a,Ns2p_a,Nop_a, yptsi_,nel_,tec_
  num_emiss = 8
  rayleighs = dblarr(num_emiss)
  epsilon_a = dblarr(num_emiss)
  triangulate,nel_,tec_,tr

  for j=0,num_emiss -1 do epsilon_a(j) = griddata(  nel_, tec_,reform(yptsi_(*,j)),xout = [nec_a], yout = [tec_a],/linear,triangles = tr )

  rayleighs(0) = 1d-6*Ns2p_a*epsilon_a(0)
  rayleighs(1) = 1d-6*Nop_a*epsilon_a(1)
  rayleighs(2) = 1d-6*Nop_a*epsilon_a(2)
  rayleighs(3) = 1d-6*Nsp_a*epsilon_a(3)
  rayleighs(4) = 1d-6*Nsp_a*epsilon_a(4)
  rayleighs(5) = 1d-6*Ns2p_a*epsilon_a(5)
  rayleighs(6) = 1d-6*Nsp_a*epsilon_a(6)
  rayleighs(7) = 1d-6*Nsp_a*epsilon_a(7)

return, rayleighs
end


function apo_deviates_cm3_using_emiss_table_fitter, p,x=x,y=y,err=err
  COMMON block1,yptsi_in_2dto1d,nel2dto1d,tec2dto1d

;;idx_want =[0,1,2,3,5,6,7] ; remove 4th of 7 indicies don't fit 4076 SII or S^+ emission line too faint at least on dusk outer edge... 
 ; model = model[*,idx_want]
;APO_cm3_model_given_emiss_table,Tec_a,nec_a,Nsp_a,Ns2p_a,Nop_a, yptsi_,nel_,tec_
 
model = APO_cm3_model_given_emiss_table(p[0],p[1],p[2],p[3],p[4], yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

idxwantz = [0,1,2,3,5,6,7]
model = model[idxwantz]

result = (y - model)/err


return,result

end




pro dawn_cm3_fitter_only_7_lines_using_tables

  restore, filename ='grids.sav', /verbose

  ;COMMON block1, i, n_in, Te_in, xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, x_grid, x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,te_noww,pp
  ;COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,te_fixed,nec_fixed_found
 ; COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed
 ; COMMON block2,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne,rhoORne,rhocne
  ;COMMON block1,x_grid,x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,Te_fixed,nec_fixed_found,rhoRne_,rhoORne_,rhocne_,Bin_,Bout_,nec_pl,totchmix_for_apo,nspmix,ns2pmix,nopmix,lambdaa ;print,min(x_grid),max(x_grid)
  COMMON block1,yptsi_in_2dto1d,nel2dto1d,tec2dto1d
 
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
  
  ;write_csv,'x_grid_full.csv',x_grid
 ;write_csv,'dawn_grid_full_unsmoothed.csv',dawn_grid
 ; write_csv,'err_dawn_grid_full_unsmoothed.csv',err_dawn_grid
 ;  write_csv,'dusk_grid_full_unsmoothed.csv',dusk_grid
 ;  write_csv,'err_dusk_grid_full_unsmoothed.csv',err_dusk_grid
    
  ;  stop
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
 ; dawn_grid = dusk_grid
  ;err_dawn_grid = err_dusk_grid
  XTITTLEE = 'Dawn $\rho_c$ ($R_J$)'
    for jj = 0, 7 do begin
    dawn_grid(*,jj) = reverse(dawn_grid(*,jj))
    err_dawn_grid(*,jj) =reverse(err_dawn_grid(*,jj))
  endfor

  ;for llll = 0, 5 do begin


    num_of_fit_params = 2
    p_output1=dblarr(n_x_grid,num_of_fit_params)
    p_error_output1=dblarr(n_x_grid,num_of_fit_params)
    yptsi_output1=dblarr(n_x_grid,8)



    ;For i= 0, n_elements(x_grid)-1 do begin

    ;x0 = x_grid(n_elements(x_grid)-1-i) ;- x_step_size/2d


    ;if ( x0 gt 5.7d) then te_noww = interpol(bag94_te_for_dusk(1,*),bag94_te_for_dusk(0,*),x0) else te_noww = 3.79556615d*(x0/5.7d)^(10.1797d) ; (for dusk)

    ; if ( x0 le 5d) then te_noww = 1d ; (for dusk)


    ;  if ( x0 gt 5.86d) then te_noww = interpol(bag94_te_for_dusk(1,*),bag94_te_for_dusk(0,*) + 0.16 ,x0) else te_noww = 3.79556615d*(x0/5.86d)^(10.485d) ;(for dawn)

    ; if ( x0 le 5.16d) then te_noww = 1d ;(for dawn)

    ; if ( x0 gt 5.9d) then te_noww = interpol(bag94_te_for_dusk(1,*),bag94_te_for_dusk(0,*) + 0.2d ,x0) else te_noww = 3.79556615d*(x0/5.9d)^(10.5614d) ;(for dawn)

    ; if ( x0 le 5.2d) then te_noww = 1d ;(for dawn)


    ;  slit_pos_vec=[x0,y0,z0]

    ;y=[dawn_grid(n_elements(x_grid)-1-i,0),dawn_grid(n_elements(x_grid)-1-i,1),dawn_grid(n_elements(x_grid)-1-i,2),dawn_grid(n_elements(x_grid)-1-i,3),dawn_grid(n_elements(x_grid)-1-i,4),dawn_grid(n_elements(x_grid)-1-i,5),dawn_grid(n_elements(x_grid)-1-i,6),dawn_grid(n_elements(x_grid)-1-i,7)]
    ;err=[err_dawn_grid(n_elements(x_grid)-1-i,0),err_dawn_grid(n_elements(x_grid)-1-i,1),err_dawn_grid(n_elements(x_grid)-1-i,2),err_dawn_grid(n_elements(x_grid)-1-i,3),err_dawn_grid(n_elements(x_grid)-1-i,4),err_dawn_grid(n_elements(x_grid)-1-i,5),err_dawn_grid(n_elements(x_grid)-1-i,6),err_dawn_grid(n_elements(x_grid)-1-i,7)]

    ;x = [3722.,3726.,3729.,4069.,4076.,6312.,6716.,6731.];
    ;x = [3726.,3729.,4069.,6312.,6716.,6731.];
    ;x = [4069.d,4076.d,6716.d,6731.d]; all 4 sp lines including low res one...
    ;x = [4069.d,6716.d,6731.d]; 3 sp lines without low res one...

    ;y=[dawn_grid(n_elements(x_grid)-1-i,1),dawn_grid(n_elements(x_grid)-1-i,2),dawn_grid(n_elements(x_grid)-1-i,3),dawn_grid(n_elements(x_grid)-1-i,5),dawn_grid(n_elements(x_grid)-1-i,6),dawn_grid(n_elements(x_grid)-1-i,7)]
    ;err=[err_dawn_grid(n_elements(x_grid)-1-i,1),err_dawn_grid(n_elements(x_grid)-1-i,2),err_dawn_grid(n_elements(x_grid)-1-i,3),err_dawn_grid(n_elements(x_grid)-1-i,5),err_dawn_grid(n_elements(x_grid)-1-i,6),err_dawn_grid(n_elements(x_grid)-1-i,7)]


    ; y=[dawn_grid(n_elements(x_grid)-1-i,3),dawn_grid(n_elements(x_grid)-1-i,4),dawn_grid(n_elements(x_grid)-1-i,6),dawn_grid(n_elements(x_grid)-1-i,7)]
    ; err=[err_dawn_grid(n_elements(x_grid)-1-i,3),err_dawn_grid(n_elements(x_grid)-1-i,4),err_dawn_grid(n_elements(x_grid)-1-i,6),err_dawn_grid(n_elements(x_grid)-1-i,7)]



    idx_want=[0,1,2,3,5,6,7];[3,6,7]
    dawn_grid_temp = dawn_grid(*,idx_want);so currently 7 lines ;
    ; err_dawn_grid(*,3) = err_dawn_grid(*,3)
    err_dawn_grid_temp = err_dawn_grid(*,idx_want);;so currently 7 lines (*,idx_want)


    ;y needs to be setup to for [model] =  reform(dblarr(num_elements,8),num_elements*8)
    ndawn_grid = n_elements(dawn_grid_temp)
   ; y=reform(dawn_grid_temp,ndawn_grid);[dawn_grid(n_elements(x_grid)-1-i,3),dawn_grid(n_elements(x_grid)-1-i,6),dawn_grid(n_elements(x_grid)-1-i,7)]
   ; err=reform(err_dawn_grid_temp,ndawn_grid);[50d*err_dawn_grid(n_elements(x_grid)-1-i,3),err_dawn_grid(n_elements(x_grid)-1-i,6),err_dawn_grid(n_elements(x_grid)-1-i,7)]

    ;write_csv,'y_simult_whole_dusk_all8lines_APO_coadded_5runningavg_4.6-7.csv',y
    ;write_csv,'err_simult_whole_dusk_all8lines_APO_coadded_5runningavg_4.6-7.csv',err

    ;write_csv,'y_simult_whole_dawn_all8lines_APO_coadded_5runningavg_4.6-7.csv',y
    ;write_csv,'err_simult_whole_dawn_all8lines_APO_coadded_5runningavg_4.6-7.csv',err
    ;stop
   ; x = replicate(1d,ndawn_grid)
    ;50d*err_dawn_grid(n_elements(x_grid)-1-i,3); artifically enhance errobars on 4069 to get better fit in 6731 and 6716, off for now


    ;Tec_a = p[0:nx_grid - 1]
    ;nec_a = p[nx_grid:2*nx_grid - 1]
    ;nsp_a = p[2*nx_grid:3*nx_grid - 1]
    ;ns2p_a = p[3*nx_grid:4*nx_grid - 1]
    ;nop_a = p[4*nx_grid:5*nx_grid - 1]


    openr,1,'p_output_first_try_0.04res_5runningavg_fit_whole_simult_8lines_all5_params.csv'
    p_output_prev=dblarr(n_x_grid,5) ; ne x tec x discrete wavelength centers of each emission line
    readf,1,p_output_prev
    close,1


    openr,1,'necsmooth1dusk.csv'
    nec_funcform_test=dblarr(60) ; ne x tec x discrete wavelength centers of each emission line
    readf,1,nec_funcform_test
    close,1

    openr,1,'nspsmooth1dusk.csv'
    nsp_funcform_test=dblarr(60) ; ne x tec x discrete wavelength centers of each emission line
    readf,1,nsp_funcform_test
    close,1

    nec_funcform_test = reverse(nec_funcform_test)
    nsp_funcform_test = reverse(nsp_funcform_test)
    
    
    restore,filename='expectation_nerneymodelv1_ceq_radial_model_nsp+2ns2p+nop_over_nec.sav',/verbose
    ;save,rho_ceq_for_totchmix ,totchmix_for_apo,filename='expectation_nerneymodelv1_ceq_radial_model_nsp+2ns2p+nop_over_nec.sav',/verbose

totchmix_for_apo = interpol(totchmix_for_apo,rho_ceq_for_totchmix + 0.24d,x_grid)

    
   tec = replicate(3.8d,n_elements(x_grid))
    
    nec = replicate(1200d,n_elements(x_grid))


    Nsp = replicate(2d12,n_elements(x_grid))
    
    Ns2p = replicate(1d13,n_elements(x_grid))


    Nop = replicate(1d13,n_elements(x_grid))
    
    
    tec_error = replicate(0d,n_elements(x_grid))

    nec_error = replicate(0d,n_elements(x_grid))


    Nsp_error = replicate(0d,n_elements(x_grid))

    Ns2p_error = replicate(0d,n_elements(x_grid))


    Nop_error = replicate(0d,n_elements(x_grid))
    

    
 
modelzz = dblarr(n_elements(x_grid),8)

for q=0, n_elements(x_grid) -1  do begin
modelzz[q,*] = APO_cm3_model_given_emiss_table(Tec[q],nec[q],Nsp[q],Ns2p[q],Nop[q], yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

endfor

x_grid_higher_res = x_grid
yptsi_output1 = modelzz
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

;stop

;APO_cm3_model_given_emiss_table,Tec_a,nec_a,Nsp_a,Ns2p_a,Nop_a, yptsi_,nel_,tec_
for q=0, n_elements(x_grid) -1 do begin



   
;stop

;APO_cm3_model_given_emiss_table(p[0],p[1],p[2],p[3],p[4], yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
 tecupval = 8d
n_x_grid = n_elements(x_grid)
nx_grid = n_x_grid
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 5)
parinfo[*].value = [tec[q],nec[q],Nsp[q],Ns2p[q],Nop[q]];
if (q>0) then begin
  parinfo[*].value = [tec[q-1],nec[q-1],Nsp[q-1],Ns2p[q-1],Nop[q-1]];
  if (tec[q-1] ge tecupval) then  parinfo[*].value = [3.8d,2200d,Nsp[q-1],Ns2p[q-1],Nop[q-1]];
endif


;if ((q>0) and (tec[q-1] = 6d)) then parinfo[*].value = [5d,nec[q-1],Nsp[q-1],Ns2p[q-1],Nop[q-1]];


parinfo[0].limited[0] = 1
parinfo[0].limits[0]  = 1d
parinfo[0].limited[1] = 1
parinfo[0].limits[1]  = tecupval


parinfo[1].limited[0] = 1
parinfo[1].limits[0]  = 100d
parinfo[1].limited[1] = 1
parinfo[1].limits[1]  = 6000d


parinfo[2].limited[0] = 1
parinfo[2].limits[0]  = 1d11
parinfo[2].limited[1] = 1
parinfo[2].limits[1]  = 1d15


parinfo[3].limited[0] = 1
parinfo[3].limits[0]  = 1d11
parinfo[3].limited[1] = 1
parinfo[3].limits[1]  = 1d15


parinfo[4].limited[0] = 1
parinfo[4].limits[0]  = 1d11
parinfo[4].limited[1] = 1
parinfo[4].limits[1]  = 1d15

idxwantzsss = [0,1,2,3,5,6,7]


y = reform(dawn_grid_temp(q,*))

err = reform(err_dawn_grid_temp(q,*))

x = replicate(1d,n_elements(y))

fa = {X:x, Y:y, ERR:err}


; if llll = 0 then stop;stop
p=mpfit('apo_deviates_cm3_using_emiss_table_fitter',functargs=fa,parinfo=parinfo,perror=perror);,ftol=1d-40,gtol=1d-40,xtol=1d-40)
print,p
print,perror

tec[q] = p[0]
nec[q] = p[1]
Nsp[q] = p[2]
Ns2p[q] = p[3]
Nop[q] = p[4]

tec_error[q] = perror[0]
nec_error[q] = perror[1]
Nsp_error[q] = perror[2]
Ns2p_error[q] = perror[3]
Nop_error[q] = perror[4]


endfor



Nec_column_given_ratio = (Nsp + 2d*Ns2p + Nop)/totchmix_for_apo


Nec_column_given_ratio_error = Sqrt(Nsp_error^2d + 4d*Ns2p_error^2d + Nop_error^2d)/totchmix_for_apo

LOS_distance_given_ratio = (Nec_column_given_ratio/nec)/(7.1492d9) ; in RJ given 1 RJ = 71492 km = 7.1492d7 m = 7.1492d9 cm
LOS_distance_given_ratio_error = ((Nec_column_given_ratio/nec)*Sqrt((Nec_column_given_ratio_error/Nec_column_given_ratio)^2d + (nec_error/nec)^2d))/(7.1492d9) 

Nspmix = Nsp/Nec_column_given_ratio 
Ns2pmix = Ns2p/Nec_column_given_ratio 
Nopmix = Nop/Nec_column_given_ratio 

Nspmix_error = Nspmix*Sqrt((Nsp_error/Nsp)^2d + (Nec_column_given_ratio_error/Nec_column_given_ratio)^2d)
Ns2pmix_error = Ns2pmix*Sqrt((Ns2p_error/Ns2p)^2d + (Nec_column_given_ratio_error/Nec_column_given_ratio)^2d)
Nopmix_error = Nopmix*Sqrt((Nop_error/Nop)^2d + (Nec_column_given_ratio_error/Nec_column_given_ratio)^2d)


modelz = dblarr(n_elements(x_grid),8)

for q=0, n_elements(x_grid) -1  do modelz[q,*] = APO_cm3_model_given_emiss_table(Tec[q],nec[q],Nsp[q],Ns2p[q],Nop[q], yptsi_in_2dto1d,nel2dto1d,tec2dto1d)



x_grid_higher_res = x_grid
yptsi_output1 = modelz
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



leg22.save,'spec_sim_cm3_fit_only_7_lines_emiss_tables_DAWN_coadded_0.04res_5pointavg.png',resolution=300

;doubles_L2_regularization_lambdaeq0.00017_let_it
;0it


left = 0.2
right = 0.05
MARGINz0 =[left,0.1,right,0.1]
MARGINz1 =[left,0.1,right,0.0]
MARGINz2 =[left,0.25,right,0.0]


p1=errorplot(x_grid,tec,tec_error,ytitle='$T_{ec}$ (eV)',layout=[2,4,1],margin =   MARGINz0,XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[1,10]);,/ylog)
p2=errorplot(x_grid,nec,nec_error,ytitle='$n_{ec}$ ($cm^{-3}$)',layout=[2,4,2],/current,margin =   MARGINz0,XTICKFORMAT="(A1)",yrange=[1000,4000],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p3=errorplot(x_grid,Nsp,Nsp_error,ytitle='$n_{S^{+}}$ ($cm^{-3}$)',layout=[2,4,3],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[1d12,2d13],/ylog,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p4=errorplot(x_grid,Nspmix,Nspmix_error,ytitle='$n_{S^{+}}/n_{ec}$ ',layout=[2,4,4],/current,margin =   MARGINz1,yrange=[0.01,0.8],XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p5=errorplot(x_grid,Ns2p,Ns2p_error,ytitle='$n_{S^{++}}$ ($cm^{-3}$)',layout=[2,4,5],/current,margin =   MARGINz1,XTICKFORMAT="(A1)",yrange=[8d12,5d13],/ylog,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p6=errorplot(x_grid,Ns2pmix,Ns2pmix_error,ytitle='$n_{S^{++}}/n_{ec}$ ',layout=[2,4,6],/current,margin =   MARGINz1,yrange=[0.01,0.5],XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p7=errorplot(x_grid,Nop,nop_error,ytitle='$n_{O^{+}}$ ($cm^{-3}$)',layout=[2,4,7],/current,margin =   MARGINz2,yrange=[8d12,5d13],/ylog,xtitle=   XTITTLEE,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)
p8=errorplot(x_grid,Nopmix,Nopmix_error,xtitle=   XTITTLEE,ytitle='$n_{O^{+}}/n_{ec}$ ',layout=[2,4,8],/current,margin =   MARGINz2,yrange=[0.01,0.9],SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01);,/ylog)

;print,min(te_fixed),min(nec_fixed_found),min(nsp_fixed_found),min(ns2p_fixed_found),min(nop_fixed_found)


p8.save,'given_fixed_sum_ratio_outs_plots_cm3_fit_only_7_lines_emiss_tables_DAWN_coadded_0.04res_5pointavg.png',resolution=300


left = 0.1
right = 0.05
MARGINz0 =[left,0.1,right,0.1]
MARGINz1 =[left,0.1,right,0.0]
MARGINz2 =[left,0.15,right,0.0]


p11112=errorplot(x_grid,Nec_column_given_ratio,Nec_column_given_ratio_error,ytitle='$N_{ec}$ ($cm^{-2}$)',layout=[1,2,1],margin =   MARGINz0,XTICKFORMAT="(A1)",SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[1d13,2d14],/ylog)
p11113=errorplot(x_grid,LOS_distance_given_ratio,LOS_distance_given_ratio_error,ytitle='LOS Distance = $N_{ec}$/$n_{ec}$  ($R_J$)',layout=[1,2,2],/current,margin =   MARGINz2,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,xtitle=   XTITTLEE)



p11113.save,'given_fixed_sum_ratio_LOS_implied_distance_cm3_fit_only_7_lines_emiss_tables_DAWN_coadded_0.04res_5pointavg.png',resolution=300


;p999 = errorplot(x_grid,totchmix,totchmix_error,ytitle='$\Sigma Z_i n_i $ /$n_{ec}$',xtitle='Dusk $\rho_c$ ($R_J$)',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,1.5],name='Fit');,/ylog,yrange=[0.01,1.5] )
 ; p9999 = plot(x_grid ,totchmix_for_apo,LINESTYLE=2,name='Nominal Model',/overplot,xrange=[4.5,7.]);,/ylog,yrange=[0.01,1.5] )

;leg = LEGEND(TARGET=[ p999, p9999], POSITION=[6.5,1.4], $
;  /DATA, /AUTO_TEXT_COLOR)
;leg.save, 'new_comp_doubles_L2_regularization_lambdaeq' + lambdavalstring + '_let_it_neutrality_check_7lines_17_params_func_form_at_once_DAWN_coadded_0.04res_5pointavg_tec_pl_nec_from_ratio_then_fixed_mixr_iterate_tec_nec_layers.png',resolution=300
;doubles_L2_regularization_lambdaeq0.000000017_let_it
;0it

model = modelz

save,tec,tec_error,nec,nec_error,Nsp,Nsp_error,Ns2p,Ns2p_error,Nop,Nop_error,Nspmix,Nspmix_error,Ns2pmix,Ns2pmix_error,Nopmix,Nopmix_error,Nec_column_given_ratio,Nec_column_given_ratio_error,LOS_distance_given_ratio,LOS_distance_given_ratio_error,x_grid,dawn_grid,err_dawn_grid,model, FILENAME = 'save_file_cm3_fit_only_7_lines_emiss_tables_DAWN_coadded_0.04res_5pointavg.sav',/verbose


 stop


end



