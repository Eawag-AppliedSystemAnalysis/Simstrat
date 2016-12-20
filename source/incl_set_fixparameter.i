! ************************************************************
!     writes values of fix parameters into common blocks
! ***********************************************************

!      *** Constants  for freshwater ***
          rho_0=1000	 		! mean density for freshwater (1023 for sea water)
          cp=4182   		! mean value for freshwater (for sea water) 
!      *** General constants
          rho_air=1.2       ! density of air
          g=9.81            ! accelaration of earch     
          kappa=0.41        ! ?
          K_s=0.05          ! ?
          z0=0.5            ! is this still used?
          pi=3.141592654
!      *** parameter for k-eps model ***
          ce1=1.44
          ce2=1.92
          ce3=-0.4	
          sig_k=1.0	
          sig_e=1.3
          cmue=0.09
!      *** parameter for MY model ***
          sl=0.2		
          e1=1.8		
          e2=1.33		
          a1=0.92		
          a2=0.74		
          B1=16.6		
          B2=10.1		
          c1=0.08	
!      *** further parameters controlling water dynamic ***
          Prndtl=0.8    ! Parndl number
          k_min=1.0e-9	! minimum value allowed for turbulent kinetic energy		(should not be changed)
          epsmin=1.0e-30	! minimum value allowed for energy dissipation
          avhmin=1.0e-8	!minimum value allowed for 

!      *** obsolete parameters ? ***
          cl=1.0		
          ModSal=1	
          SalRel=172800		
