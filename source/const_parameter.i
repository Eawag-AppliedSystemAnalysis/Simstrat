      integer mxl,noout
      parameter(mxl=2000)
      parameter(noout=14)

      double precision depth,rho_0,cp,kappa,K_s
      double precision cm0,z0,Cori,ce1,ce2,ce3,dt,t_start,t_end
	  double precision pi,k_min,epsmin,avhmin,g
      double precision sig_e,Prndtl,CD,C10,cde,cmue
      double precision Sl,e1,e2,a1,a2,b1,b2,c1,cl,SalRel
      double precision alpha_Seiche, CDeff, q_NN, rho_air,sig_k
      double precision H_A, H_K, H_V, H_W
      double precision p_radin,p_windfun
      double precision zs(0:1000)	 

      double precision L_min

      double precision fgeo,fgeo_add(0:mxl),fsed,fsed_add(0:mxl)
      double precision air_press

      integer disp_diagnose,disp_simulation


      integer stab,Mod,fluxcond,NBC,xl,xl_ini,num_save,depth_save,modsal,pgrad
      integer igoal,igoal_s(1:noout)
      integer xlk,xlu,indexk_save(mxl),indexu_save(mxl) 
      integer save_vctr,ifit, norm_ctr
      character*100 ParName,InitName, MorphName,ForcingName
      character*100 AbsorpName,OutName,GridName,zoutName,toutName
      character*100 QinpName,QoutName,TinpName,SinpName
	  integer adv

      integer salctr,delsal 	 
   
      common /com_depth/  depth
      common /com_pi/     pi
      common /com_rho_0/  rho_0
      common /com_rho_air/rho_air
      common /com_cp/     cp
      common /com_g/      g
      common /com_xl/     xl,xl_ini
      common /com_cde/    cde
      common /com_cmue/   cmue
      common /com_Cori/   Cori
      common /com_dt/     dt
      common /com_sig_e/  sig_e
      common /com_sig_k/  sig_k
      common /com_Prndtl/ Prndtl
      common /com_z0/     z0
      common /com_k_min/  k_min,L_min
      common /com_epsmin/ epsmin
      common /com_avhmin/ avhmin
      common /com_NBC/    NBC
      common /com_pgrad/  pgrad
      common /com_stab/   stab
      common /com_Mod/    Mod
      common /com_modsal/ modsal
      common /com_fluxcond/ fluxcond
      common /com_CD/     CD,C10
      common /com_Sl/     Sl
      common /com_cm0/    cm0
      common /com_SalRel/ SalRel
      common /com_time/   t_start,t_end
      common /com_save/   num_save,depth_save
      common /com_Seiche/ alpha_Seiche,CDeff,q_NN
      common /com_dissip/ ce1,ce2,ce3
      common /com_K_s/    K_s
      common /com_kappa/  kappa
      common /com_e/      e1,e2
      common /constants/  a1,a2,b2,c1
      common /com_cl/     cl
      common /com_b1/     b1
      common /com_heat_budget/ H_A, H_K, H_V, H_W
      common /filenames/  ParName,InitName,ForcingName,MorphName,AbsorpName,OutName,GridName,
     &	  zoutName,toutName,QinpName,QoutName,TinpName,SinpName, adv
      common /savectr/    xlk,xlu,indexk_save,indexu_save
      common /com_igoal/  igoal,igoal_s
      common /savez/      zs,save_vctr
      common /heat_fit/   p_radin,p_windfun
      common /geo_therm/ fgeo,fgeo_add,fsed,fsed_add
      common /air_p/ air_press
      common /c_ctr/ norm_ctr

      common /buoy/       salctr,delsal  
	  
      common /display/        disp_diagnose,disp_simulation

      common /fit_controle/  ifit
