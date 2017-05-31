module SimstratOutput
    use SimstratModel
    use utilities
    implicit none

    private
    public save_ini, write_out, prep_outav, write_out_new, write_text, avstate

contains
    !Prepare parameters for output of the results
    !####################################################################
    subroutine save_ini()
    !####################################################################

        implicit none


        real(dp) :: t_out(0:9000), test(0:nz_max)
        integer          :: i,j,tctr

        ! Read z-vector for vertical output
        if(disp_dgn/=0) write(6,*) 'Output depths  : ',trim(zoutName)
        open(15,status='old',file=zoutName)
        read(15,*)
        do i=0,nz_max
            read(15,*,end=8) zsave(i)
        end do
    8   if(i==nz_max) write(6,*) 'Only first ',nz_max,' values of file read.'
        close(15)
	write(*,*) zsave
        nsave=i-1

        if(OutBin) write(80) igoal
        if(OutBin) write(80) nsave

        if (nsave==0) then             ! Take every nth index for output
            j=0
            depth_save=int(zsave(0))
            do i=nz,0,-depth_save
                j=j+1
                index_upp_save(j)=i
                index_cent_save(j)=i
            end do
            nz_upp=j
            nz_cent=j

            if (index_upp_save(nz_upp)>0) then
                nz_upp=nz_upp+1
                index_upp_save(nz_upp)=0
            end if
            if (index_cent_save(nz_cent)>1) then
                nz_cent=nz_cent+1
                index_cent_save(nz_cent)=1
            end if

            ! Write vertical vector to file
            if (OutBin) then
                if (abs(igoal)>9) then
                    write(80) nz_cent
                    write(80) (z_upp(nz)-z_cent(index_cent_save(i)),i=1,nz_cent)
                    write(80) nz_upp
                    write(80) (z_upp(nz)-z_upp(index_upp_save(i)),i=1,nz_upp)
                elseif (abs(igoal)<5) then
                    write(80) nz_cent
                    write(80) (z_upp(nz)-z_cent(index_cent_save(i)),i=1,nz_cent)
                elseif (abs(igoal)<=8) then
                    write(80) nz_upp
                    write(80) (z_upp(nz)-z_upp(index_upp_save(i)),i=1,nz_upp)
                else
                    write(6,*) 'Error: wrong choice of index to goal parameter.'
                    stop
                end if
            end if
        else                               ! If depth vector is given by user
            if (OutBin) then
                write(80) nsave+1
                write(80) (zsave(i),i=0,nsave)
            end if
            !Output depths are absolute points
            test(0:nsave) = zsave(0:nsave)
            do i=0,nsave
                zsave(nsave-i) = z_zero + test(i)
            end do
	    write(*,*) "------------------------------__"
	    write(*,*) zsave
            
        end if
        if (disp_dgn/=0) then
            write(6,*) 'Successfully read'
            write(6,*)
        end if

        ! Read t-vector for temporal output
        if(disp_dgn/=0) write(6,*) 'Output times   : ',trim(toutName)
        open(15,status='old',file=toutName)
        read(15,*)
        do i=0,9000
            read(15,*,end=9) t_out(i)
        end do
    9   if(i==9000) write(6,*) 'Only first 9000 values of file read'
        close(15)
        tctr=i-1

        if (tctr/=0) then                  ! Indices to save
            if (t_start>t_out(0)) then
                write(6,*) 'Error: simulation start time is larger than first output time'
                stop
            end if

            tout_ctr1(0)=(t_out(0)-t_start)*86400/dt
            if (int(tout_ctr1(0))==0) then
                tout_ctr2(0)=(t_out(0)-t_start)*86400
                test(0)=t_start+tout_ctr2(0)
            else
                tout_ctr2(0)=((t_out(0)-t_start)*86400)/int(tout_ctr1(0))
                test(0)=t_start
                do j=1,int(tout_ctr1(0))
                    test(0)=test(0)+tout_ctr2(0)/86400
                end do
            end if

            do i=1,tctr
                tout_ctr1(i)=(t_out(i)-test(i-1))*86400/dt
                if (int(tout_ctr1(i))==0) then
                    tout_ctr2(i)=(t_out(i)-test(i-1))*86400
                    test(i)=test(i-1)+tout_ctr2(i)
                    write(6,*) 'Warning: time interval for output is smaller than dt for iteration'
                else
                    tout_ctr2(i)=((t_out(i)-test(i-1))*86400)/int(tout_ctr1(i))
                    test(i)=test(i-1)
                    do j=1,int(tout_ctr1(i))
                        test(i)=test(i)+tout_ctr2(i)/86400
                    end do
                end if
            end do
            t_end=test(tctr)
            write_tout=0
        else                               ! Interval given for output
            write_tout=int(t_out(0))
        end if
        if (disp_dgn/=0) then
            if(write_tout/=0) write(6,*) 'Interval [days]: ',write_tout*dt/86400
            write(6,*) 'Successfully read'
            write(6,*)
        end if

        return
    end subroutine save_ini

    !####################################################################
    !subroutine write_out(datum,u,v,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,z_cent,M)
    subroutine write_out(datum,U,V,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche)
    !####################################################################
        implicit none

        ! Global Declarations
        real(dp), intent(in) :: datum
        real(dp), intent(in) :: U(0:nz),V(0:nz),T(0:nz),S(0:nz)
        real(dp), intent(in) :: k(0:nz),eps(0:nz),num(0:nz),nuh(0:nz)
        real(dp), intent(in) :: B(0:nz),P(0:nz),NN(0:nz)
        real(dp), intent(in) :: P_Seiche(0:nz),E_Seiche
        !real(dp) M(0:nz_max)

        ! Local Declarations
        real(dp) :: xs(0:nz_max)
        integer :: i!,ist

        write(80) datum
        if (nsave==0) then
            if (igoal_s(1)==1) then
                write(80) (U(index_cent_save(i)),i=1,nz_cent)
            end if
            if (igoal_s(2)==1) then
                write(80) (V(index_cent_save(i)),i=1,nz_cent)
            end if
            if (igoal_s(3)==1) then
                write(80) (T(index_cent_save(i)),i=1,nz_cent)
            end if
            if (igoal_s(4)==1) then
                write(80) (S(index_cent_save(i)),i=1,nz_cent)
            end if
            if (igoal_s(5)==1) then
                write(80) (k(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(6)==1) then
                write(80) (eps(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(7)==1) then
                write(80) (num(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(8)==1) then
                write(80) (nuh(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(9)==1) then
                write(80) (B(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(10)==1) then
                write(80) (P(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(11)==1) then
                write(80) (P_Seiche(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(12)==1) then
                write(80) (NN(index_cent_save(i)),i=1,nz_upp)
            end if
            if (igoal_s(13)==1) then
                write(80) E_Seiche
            end if
        else
            if (igoal_s(1)==1) then
                call Interp(z_cent,U,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(2)==1) then
                call Interp(z_cent,V,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(3)==1) then
                call Interp(z_cent,T,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(4)==1) then
                call Interp(z_cent,S,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(5)==1) then
                call Interp(z_cent,k,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(6)==1) then
                call Interp(z_cent,eps,nz,zsave,xs,nsave)
                 write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(7)==1) then
                call Interp(z_cent,num,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(8)==1) then
                call Interp(z_cent,nuh,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(9)==1) then
                 call Interp(z_cent,B,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(10)==1) then
                call Interp(z_cent,P,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(11)==1) then
                call Interp(z_cent,P_Seiche,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(12)==1) then
                call Interp(z_cent,NN,nz,zsave,xs,nsave)
                write(80) (xs(i),i=nsave, 0, -1)
            end if
            if (igoal_s(13)==1) then
                write(80) E_Seiche
            end if
        endif     ! save -ctr

        return
    end

    !####################################################################
    subroutine prep_outav(i_av,U_av,V_av,T_av,S_av,k_av,eps_av,num_av,nuh_av,B_av,P_av,NN_av,P_Seiche_av,E_Seiche_av)
    !####################################################################
        implicit none

        ! Global Declarations
        real(dp), intent(inout) :: U_av(0:nz),V_av(0:nz),T_av(0:nz),S_av(0:nz)
        real(dp), intent(inout) :: k_av(0:nz),eps_av(0:nz),num_av(0:nz),nuh_av(0:nz)
        real(dp), intent(inout) :: B_av(0:nz),P_av(0:nz),NN_av(0:nz)
        real(dp), intent(inout) :: P_Seiche_av(0:nz),E_Seiche_av
        integer, intent(in) :: i_av

        if(igoal_s(1)==1) U_av(0:nz)=U_av(0:nz)/i_av
        if(igoal_s(2)==1) V_av(0:nz)=V_av(0:nz)/i_av
        if(igoal_s(3)==1) T_av(0:nz)=T_av(0:nz)/i_av
        if(igoal_s(4)==1) S_av(0:nz)=S_av(0:nz)/i_av
        if(igoal_s(5)==1) k_av(0:nz)=k_av(0:nz)/i_av
        if(igoal_s(6)==1) eps_av(0:nz)=eps_av(0:nz)/i_av
        if(igoal_s(7)==1) num_av(0:nz)=num_av(0:nz)/i_av
        if(igoal_s(8)==1) nuh_av(0:nz)=nuh_av(0:nz)/i_av
        if(igoal_s(9)==1) B_av(0:nz)=B_av(0:nz)/i_av
        if(igoal_s(10)==1) P_av(0:nz)=P_av(0:nz)/i_av
        if(igoal_s(11)==1) P_Seiche_av(0:nz)=P_Seiche_av(0:nz)/i_av
        if(igoal_s(12)==1) NN_av(0:nz)=NN_av(0:nz)/i_av
        if(igoal_s(13)==1) E_Seiche_av=E_Seiche_av/i_av

        return
    end


    !####################################################################
    subroutine write_out_new(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,Q_vert)
    !####################################################################
        implicit none

        ! Global Declarations
        real(dp), intent(in) :: datum
        real(dp), intent(in) :: U(0:nz),V(0:nz),T(0:nz),S(0:nz)
        real(dp), intent(in) :: k(0:nz),eps(0:nz),nuh(0:nz)
        real(dp), intent(in) :: B(0:nz),P(0:nz),NN(0:nz)
        real(dp), intent(in) :: P_Seiche(0:nz),E_Seiche
        real(dp), intent(in) :: Q_vert(0:nz)

        ! Local Declarations
        real(dp) :: xs(0:nz_max)
        real(dp) :: zext(0:nz)
        integer :: i

        !External depths: bottom, layer centers, surface
        zext(0) = z_upp(0)
        zext(1:nz-2) = z_cent(2:nz-1)
        zext(nz-1) = z_upp(nz)

        write(80) datum

        call Interp_nan(zext,U,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(zext,V,nz-1,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(zext,T,nz-1,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(zext,S,nz-1,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)

        call Interp_nan(zext,Q_vert,nz-1,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)

        call Interp_nan(z_upp,k,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call   Interp_nan(z_upp,eps,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(z_upp,nuh,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(z_upp,B,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(z_upp,P,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(z_upp,P_Seiche,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)
        call Interp_nan(z_upp,NN,nz,zsave,xs,nsave)
        write(80) (xs(i),i=nsave, 0, -1)

        !Write surface value
        write(80) E_Seiche,z_upp(nz),z_cent(nz),U(nz),V(nz),T(nz),S(nz),k(nz),&
                  eps(nz),nuh(nz),B(nz),P(nz),P_Seiche(nz),NN(nz),Q_vert(nz)

        return
    end


    !Write output text files for physical variables
    subroutine write_text(datum,U,V,T,S,k,eps,nuh,B,P,NN,P_Seiche,E_Seiche,Q_vert)
        implicit none

        ! Global Declarations
        real(dp), intent(in) :: datum
        real(dp), intent(in) :: U(0:nz_grid),V(0:nz_grid),T(0:nz_grid),S(0:nz_grid)
        real(dp), intent(in) :: k(0:nz_grid),eps(0:nz_grid),nuh(0:nz_grid)
        real(dp), intent(in) :: B(0:nz_grid),P(0:nz_grid),NN(0:nz_grid)
        real(dp), intent(in) :: P_Seiche(0:nz_grid),E_Seiche
        real(dp), intent(in) :: Q_vert(0:nz_grid)

        ! Local Declarations
        real(dp) :: zext(0:nz-1)

        !External depths: bottom, layer centers, surface
        zext(0) = z_upp(0)
        zext(1:nz-2) = z_cent(2:nz-1)
        zext(nz-1) = z_upp(nz)
        call write_text_var(datum,zext,U,81)
        call write_text_var(datum,zext,V,82)
        call write_text_var(datum,zext,T,83)
        call write_text_var(datum,zext,S,84)
        call write_text_var(datum,zext,k,85)
        call write_text_var(datum,zext,eps,86)
        call write_text_var(datum,zext,nuh,87)
        call write_text_var(datum,zext,B,88)
        call write_text_var(datum,zext,P,89)
        call write_text_var(datum,zext,NN,90)
        call write_text_var(datum,zext,Q_vert,91)
        call write_text_var(datum,zext,P_Seiche,92)
        write(93,*)
        write(93,'(F10.4,$)') datum
        write(93,'(ES12.4,$)') E_Seiche

        return
    end

    !####################################################################
    subroutine write_text_var(datum,zext,var,fid)
    !####################################################################
        implicit none

        real(dp), intent(in) :: datum,zext(0:nz-1),var(0:nz_grid)

        real(dp) :: xs(0:nsave+1)
        integer :: i,fid

        call Interp_nan(zext,var,nz-1,zsave,xs,nsave)
        xs(nsave+1) = var(nz) !Surface value

        write(fid,*)
        write(fid,'(F10.4,$)') datum
        do i=0,nsave+1
            write(fid,'(ES12.4,$)') xs(i)
        end do

        return
    end subroutine

    !####################################################################
    subroutine avstate(i_av,U,V,T,S,k,eps,num,nuh,B,P,NN,P_Seiche,E_Seiche,U_av,V_av,&
                       T_av,S_av,k_av,eps_av,num_av,nuh_av,B_av,P_av,NN_av,P_Seiche_av,E_Seiche_av)
    !####################################################################
        implicit none

        ! Global Declarations
        real(dp), intent(in) :: U(0:nz),V(0:nz),T(0:nz),S(0:nz)
        real(dp), intent(in) :: k(0:nz),eps(0:nz),num(0:nz),nuh(0:nz)
        real(dp), intent(in) :: B(0:nz),P(0:nz),NN(0:nz)
        real(dp), intent(in) :: P_Seiche(0:nz),E_Seiche

        real(dp), intent(inout) :: U_av(0:nz),V_av(0:nz),T_av(0:nz),S_av(0:nz)
        real(dp), intent(inout) :: k_av(0:nz),eps_av(0:nz),num_av(0:nz),nuh_av(0:nz)
        real(dp), intent(inout) :: B_av(0:nz),P_av(0:nz),NN_av(0:nz)
        real(dp), intent(inout) :: P_Seiche_av(0:nz),E_Seiche_av
        integer, intent(in) :: i_av

        if (igoal_s(1)==1) then
            if (i_av==1) then
                U_av(0:nz) = U(0:nz)
            else
                U_av(0:nz) = U_av(0:nz) + U(0:nz)
            end if
        end if
        if (igoal_s(2)==1) then
            if (i_av==1) then
                V_av(0:nz) = V(0:nz)
            else
                V_av(0:nz) = V_av(0:nz) + V(0:nz)
            end if
        end if
        if (igoal_s(3)==1) then
            if (i_av==1) then
                T_av(0:nz) = T(0:nz)
            else
                T_av(0:nz) = T_av(0:nz) + T(0:nz)
            end if
        end if
        if (igoal_s(4)==1) then
            if (i_av==1) then
                S_av(0:nz) = S(0:nz)
            else
                S_av(0:nz) = S_av(0:nz) + S(0:nz)
            end if
        end if
        if (igoal_s(5)==1) then
            if (i_av==1) then
                k_av(0:nz) = k(0:nz)
            else
                k_av(0:nz) = k_av(0:nz) + k(0:nz)
            end if
        end if
        if (igoal_s(6)==1) then
            if (i_av==1) then
                eps_av(0:nz) = eps(0:nz)
            else
                eps_av(0:nz) = eps_av(0:nz) + eps(0:nz)
            end if
        end if
        if (igoal_s(7)==1) then
            if (i_av==1) then
                num_av(0:nz) = num(0:nz)
            else
                num_av(0:nz) = num_av(0:nz) + num(0:nz)
            end if
        end if
        if (igoal_s(8)==1) then
            if (i_av==1) then
                nuh_av(0:nz) = nuh(0:nz)
            else
                nuh_av(0:nz) = nuh_av(0:nz) + nuh(0:nz)
            end if
        end if
        if (igoal_s(9)==1) then
            if (i_av==1) then
                B_av(0:nz) = B(0:nz)
            else
                B_av(0:nz) = B_av(0:nz) + B(0:nz)
            end if
        end if
        if (igoal_s(10)==1) then
            if (i_av==1) then
                P_av(0:nz) = P(0:nz)
            else
                P_av(0:nz) = P_av(0:nz) + P(0:nz)
            end if
        end if
        if (igoal_s(11)==1) then
            if (i_av==1) then
                P_Seiche_av(0:nz) = P_Seiche(0:nz)
            else
                P_Seiche_av(0:nz) = P_Seiche_av(0:nz) + P_Seiche(0:nz)
            end if
        end if
        if (igoal_s(12)==1) then
            if (i_av==1) then
                NN_av(0:nz) = NN(0:nz)
            else
                NN_av(0:nz) = NN_av(0:nz) + NN(0:nz)
            end if
        end if
        if (igoal_s(13)==1) then
            if (i_av==1) then
                E_Seiche_av = E_Seiche
            else
                E_Seiche_av = E_Seiche_av + E_Seiche
            end if
        end if

        return
    end subroutine
end module SimstratOutput
