      integer mxl
      parameter(mxl=1000)

      double precision depth,dt,t_start,t_end
      double precision pi,g,cp,rho_0,rho_air,kappa,K_s,z0,Prndtl
      double precision k_min,eps_min,avh_min
      double precision sig_e,sig_k,ce1,ce2,ce3,cde,cmue,cm0
      double precision a1,a2,b1,b2,c1,e1,e2,sl
      double precision p_air,a_seiche,q_NN,CD,C10,Cori
      double precision p_radin,p_windf,beta_sol,albsw,f_wind
      double precision fgeo,fgeo_add(0:mxl),fsed,fsed_add(0:mxl)
      integer salctr,delsal,adv
      character*100 ParName,MorphName,InitName,ForcingName,AbsorpName
      character*100 GridName,zoutName,toutName,PathOut
      character*100 QinpName,QoutName,TinpName,SinpName
      integer Mod,Stab,ModFlux,NBC,WindFilt,ModSNorm,ModC10,ModInflow,Pgrad,ModSal
      logical OutBin
      integer xli,xl,num_save,depth_save,nsave
      double precision zsave(0:mxl)
      integer igoal,igoal_s(1:14),ifit
      integer xlk,xlu,indexk_save(1:mxl),indexu_save(1:mxl)
      integer disp_dgn,disp_sim

      common /vgrid/      depth,xli,xl
      common /time/       dt,t_start,t_end
      common /constants/  pi,g,cp,rho_0,rho_air,kappa,K_s,z0,Prndtl
      common /keps_min/   k_min,eps_min,avh_min
      common /keps_sigma/ sig_e,sig_k
      common /keps_dissip/ce1,ce2,ce3
      common /keps_cme/   cde,cmue,cm0
      common /MY_cst/     a1,a2,b1,b2,c1,e1,e2,sl
      common /settings/   Mod,Stab,ModFlux,NBC,WindFilt,ModSNorm,ModC10,ModInflow,Pgrad,ModSal
      common /phys_par/   p_air,C10,CD
      common /fit_par/    p_radin,p_windf,f_wind,beta_sol,albsw
      common /seiche_par/ a_seiche,q_NN
      common /calc_par/   Cori,adv
      common /geo_sed/    fgeo,fgeo_add,fsed,fsed_add
      common /buoy/       salctr,delsal
      common /input/      ParName,MorphName,GridName,zoutName,toutName,InitName,ForcingName,&
                          AbsorpName,QinpName,QoutName,TinpName,SinpName
      common /output/     PathOut,OutBin
      common /avg/        igoal,igoal_s,ifit
      common /savectr/    xlk,xlu,indexk_save,indexu_save
      common /saving/     num_save,depth_save
      common /savez/      zsave,nsave
      common /display/    disp_dgn,disp_sim
