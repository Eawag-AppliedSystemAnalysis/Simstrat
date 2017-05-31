    module Stability
    use SimstratModel
    implicit none


    private
    public doBuoyancy, doStabilityFunctions

contains

    !####################################################################
    subroutine doStabilityFunctions(k,eps,T,S,rho,NN,cmue1,cmue2)
    !####################################################################
        implicit none

        ! Global variables
        real(dp), intent(in) :: T(0:),S(0:)
        real(dp), intent(in) :: k(0:),eps(0:)
        real(dp), intent(inout) :: cmue1(0:),cmue2(0:),NN(0:),rho(0:)
        real(dp) :: beta_t(0:nz)
        ! Local variables
        real(dp) :: beta
        integer :: i


        call doBuoyancy(T,S,rho,NN)

        !In the domain
        do i=1,nz-1
            if (stab==1) then
                call cmue_cn(cmue1(i),cmue2(i))
            else if (stab==2) then

                beta = NN(i)*(k(i)/eps(i))**2
                beta_t(i)  = beta
                call cmue_qe(beta,cmue1(i),cmue2(i))
            end if
        end do
        beta_t(0) = 0
        beta_t(nz) = 0

        !At the boundaries
        cmue1(0) = cmue1(1)
        cmue2(0) = cmue2(1)
        cmue1(nz) = cmue1(nz-1)
        cmue2(nz) = cmue2(nz-1)

        return
    end

    !####################################################################
    subroutine doBuoyancy(T,S,rho,NN)
    !####################################################################
        !use SimstratModel, only: nz
        implicit none

        ! Global variables
        real(dp), intent(in) :: T(0:),S(0:)
        real(dp), intent(inout) :: NN(0:),rho(0:)

        ! Local variables
        real(dp) :: a(0:nz),buoy(0:nz)
        real(dp) :: rho0t(0:nz),rho0st(0:nz)
        integer :: i


        if (salctr==0) then                ! salinity is zero everywhere
            do i=1,nz-1
                a(i)= -68.0_dp+T(i)*(18.2091_dp+T(i)*(-0.30866_dp+T(i)*&
                     (5.3445e-3_dp+T(i)*(-6.0721e-5_dp+T(i)*(3.1441e-7_dp)))))
                !if (press/=0) then   ! ignore this pressure thing for alpha in first approximation
                !    a(i)= a(i) + (0.3682+T(i)*(-1.520e-2+T(i)*(1.91e-4)))*p(i)
                !end if
                a(i)= 1.0e-6_dp*a(i)
                NN(i)= g*a(i)*(meanint(i)*(T(i)-T(i+1))+(T(i)+273.15_dp)/cp)
            end do
        else
            if (delsal==0) then        ! salinity gradient is zero everywhere
                do i=1,nz-1
                    a(i)= -68.0_dp+T(i)*(18.2091_dp+T(i)*(-0.30866_dp+T(i)*&
                         (5.3445e-3_dp+T(i)*(-6.0721e-5_dp+T(i)*(3.1441e-7_dp)))))
                    a(i)= a(i) + (4.599_dp+T(i)*(-0.1999_dp+T(i)*(2.79e-3_dp)))*S(i)
                    !if (press/=0) then   ! ignore this pressure thing for alpha in first approximation
                    !    a(i) = a(i) + (0.3682+T(i)*(-1.520e-2+T(i)*(1.91e-4))-S(i)*(4.613e-3))*p(i)
                    !end if
                    a(i)= 1.0e-6_dp*a(i)
                    NN(i)= g*a(i)*(meanint(i)*(T(i)-T(i+1))+(T(i)+273.15_dp)/cp)
                end do
            else
                !if (press/=0) then
                !    rho0t= (0.9998395+t.*(6.7914e-5 +t.*(-9.0894e-6+t.*(1.0171e-7+t.*(-1.2846e-9 +t.*(1.1592e-11
                !           +t.*(-5.0125e-14)))))));
                !    kbart= (19652.17 +t.*(148.113 +t.*(-2.293 +t.*(1.256e-2+t.*(-4.18e-5)))));
                !    kbarpt= ((3.2726 +t.*(-2.147e-4 +t.*(1.128e-4))).*p);
                !    kbar= kbart+kbarpt
                !    rho0st= ((8.181e-4 +t.*(-3.85e-6 +t.*(4.96e-8))).*s);
                !    kbarspt= ((53.238 -0.313.*t +5.728e-3.*p).*s);
                !    kbar= kbar+kbarspt;
                !    rho= rho_0*(rho0t + rho0st)./(1.0 - p./kbar);
                !    if ~isempty(f)
                !        rhoS0= rho0t./(1.0 - p./(kbart+kbarpt));
                !        rho= rhoS0.*(1-fc)+fc.*rho;
                !    end
                !    rho= rho_0*rho0t./(1.0 - p./kbar);
                !else
                do i=0,nz
                    rho0t(i)= 0.9998395_dp+T(i)*(6.7914e-5_dp+T(i)*(-9.0894e-6_dp+T(i)*&
                        (1.0171e-7_dp+T(i)*(-1.2846e-9_dp+T(i)*(1.1592e-11_dp+T(i)*(-5.0125e-14_dp))))))
                    rho0st(i)= (8.181e-4_dp+T(i)*(-3.85e-6_dp+T(i)*(4.96e-8_dp)))*S(i)
                    rho(i)= rho_0*(rho0t(i)+rho0st(i))
                    !if (fc/=0) then
                    !    rho(i) = rho0t(i)*(1-fc)+fc*rho(i)
                    !end
                    buoy(i)= -g*(rho(i)-rho_0)/rho_0
                end do

                NN(1:nz-1) = meanint(1:nz-1)*(buoy(2:nz)-buoy(1:nz-1))

            end if
        end if
        NN(0)= NN(1)
        NN(nz)= NN(nz-1)

      !  write(*,*) "BUOY="
      !  write(*,*) (buoy(2:nz)-buoy(1:nz-1))
      !  write(*,*) "meanint="
      !  write(*,*) meanint(1:nz-1)
      !  write(*,*) ""
        return
    end

    ! Calculation of cmue
    !####################################################################
    subroutine cmue_cn(cmue1,cmue2)
    !####################################################################
        implicit none

        ! Global variables
        real(dp), intent(inout) :: cmue1,cmue2

        !Standard version of k-eps model
        cmue1 = cmue
        cmue2 = cmue/Prndtl
        !Burchard Version
        !cmue1=cde*cm0
        !cmue2=cde*cm0/Prndtl

        return
    end

    ! Calculation of cmue
    !####################################################################
    subroutine cmue_qe(beta,cmue1,cmue2)
    !####################################################################
        implicit none

        ! Global variables
        real(dp), intent(inout) :: cmue1,cmue2
        real(dp), intent(in) :: beta

        ! Local variables
        real(dp) :: gh,sm,sh


        gh = -cde**2*0.5_dp*beta
        if(gh> 0.02) gh = gh-(gh-0.02_dp)**2/(gh+0.0233_dp-2*0.02_dp)
        if(gh<-0.28) gh = -0.28_dp

        sm = 1.0_dp-3*c1-6*a1/b1-3*a2*gh*((b2-3*a2)*(1.0_dp-6*a1/b1)-3*c1*(b2+6*a1))
        sm = a1*sm/((1.0_dp-3*a2*gh*(6*a1+b2))*(1.0_dp-9*a1*a2*gh))
        sh = a2*(1.0_dp-6*a1/b1)/(1.0_dp-3*a2*gh*(6*a1+b2))

        cmue1 = sqrt(2.0_dp)*cde*sm
        cmue2 = sqrt(2.0_dp)*cde*sh

        return
    end

end module Stability
