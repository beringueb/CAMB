if (CTrans%NumSources>2 .and. State%CP%want_cl_2D_array) then
                        do w_ix=1,3 + State%num_redshiftwindows + num_cmb_freq*2
                            Delta1 = CTrans%Delta_p_l_k(w_ix,j,q_ix)
                            !if (w_ix>nscatter*2+1.and. w_ix2==3) then
                                !if (w_ix>3) then
                                 !   associate (Win => State%Redshift_w(w_ix - 3))
                                  !      write(*,*) 'windkind : ', Win%kind
                                   !     if (Win%kind == window_lensing) &
                                    !        Delta1 = Delta1 / 2 * ell * (ell + 1)
                                     !       write(*,*) 'HERE 0?', window_lensing
                                      !  if (Win%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                                            !want delta f/f - 2kappa;
                                            ! grad^2 = -l(l+1);
                                       !      write(*,*) 'HERE ?', Win%mag_index
                                        !    Delta1 = Delta1 + ell * (ell + 1) * &
                                         !       CTrans%Delta_p_l_k(3 + Win%mag_index + &
                                          !      State%num_redshiftwindows, j, q_ix)
                             !           end if
                             !       end associate
                             !   end if
                            !end if
                            do w_ix2=w_ix,3 + State%num_redshiftwindows + num_cmb_freq*2
                                Delta2 = CTrans%Delta_p_l_k(w_ix2, j, q_ix)
                                if (w_ix2>nscatter*2+1 .and. w_ix==3) then
                                    if (w_ix>nscatter*2+1.and. w_ix2==3) then
                                !if (w_ix2>nscatter*2+1 .or. w_ix2==3) then
                                !if (w_ix2>= 3.and. w_ix>=3) then
                                    !Skip if the auto or cross-correlation is included in direct Limber result
                                    !Otherwise we need to include the sources e.g. to get counts-Temperature correct
                                        if (CTrans%limber_l_min(w_ix2)/= 0 .and. j>=CTrans%limber_l_min(w_ix2) &
                                            .and. CTrans%limber_l_min(w_ix)/= 0 .and. j>=CTrans%limber_l_min(w_ix)) cycle
                                    end if
                                end if
                                Delta2 = CTrans%Delta_p_l_k(w_ix2, j, q_ix)
                                !if (w_ix2 > 3) then
                                !    associate (Win => State%Redshift_w(w_ix2 - 3))
                                !        if (Win%kind == window_lensing) &
                                !            Delta2 = Delta2 / 2 * ell * (ell + 1)
                                !        if (Win%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                                            !want delta f/f - 2kappa;
                                            ! grad^2 = -l(l+1);
                                !            Delta2 = Delta2 + ell * (ell + 1) * &
                                !                CTrans%Delta_p_l_k(3 + Win%mag_index + &
                                !                State%num_redshiftwindows, j, q_ix)
                                !        end if
                                !    end associate
                                !end if
                                iCl_Array(j,w_ix,w_ix2) = iCl_Array(j,w_ix,w_ix2)+Delta1*Delta2*apowers*dlnk
                            end do
                        end do
                        if (CP%CustomSources%num_custom_sources >0) then
                            do w_ix=1,3 + State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                                if (w_ix > 3 + State%num_redshiftwindows) then
                                    Delta1= CTrans%Delta_p_l_k(w_ix+State%num_extra_redshiftwindows,j,q_ix)
                                else
                                    Delta1= CTrans%Delta_p_l_k(w_ix,j,q_ix)
                                end if
                                do w_ix2=max(w_ix,3 + State%num_redshiftwindows +1), &
                                    3 + State%num_redshiftwindows +CP%CustomSources%num_custom_sources
                                    Delta2=  CTrans%Delta_p_l_k(w_ix2+State%num_extra_redshiftwindows,j,q_ix)
                                    iCl_Array(j,w_ix,w_ix2) = iCl_Array(j,w_ix,w_ix2) &
                                        +Delta1*Delta2*apowers*dlnk
                                end do
                            end do
                        end if
                    end if

                    if (CTrans%NumSources>2 ) then
                        if (CP%SourceTerms%limber_phi_lmin==0 .or.  &
                            CTrans%limber_l_min(3)== 0 .or. j<CTrans%limber_l_min(3)) then
                            iCl_scalar(j,C_Phi) = iCl_scalar(j,C_Phi) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)**2*dlnk
                            iCl_scalar(j,C_PhiTemp) = iCl_scalar(j,C_PhiTemp) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
                            iCl_scalar(j,C_PhiE) = iCl_scalar(j,C_PhiE) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
                        end if
                    end if
                end if
            end do

        end if !limber (j<= max_bessels_l_index)

        !Output l(l+1)C_l/OutputDenominator
        ctnorm=(ell*ell-1)*(ell+2)*ell
        dbletmp=(ell*(ell+1))/OutputDenominator*const_fourpi
        if (State%CP%want_cl_2D_array) then
            fac=1
            fac(2) = sqrt(ctnorm)
            if (CTrans%NumSources > 2) then
                fac(3) = sqrt(ell*(ell+1)*CP%ALens) !Changed Dec18 for consistency
                do w_ix=3 + State%num_redshiftwindows + 1,3 + State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                    !nscal= CP%CustomSources%custom_source_ell_scales(w_ix - State%num_redshiftwindows -3)
                    do i=1, nscal
                        fac(w_ix) = fac(w_ix)*(ell+i)*(ell-i+1)
                    end do
                    fac(w_ix) = sqrt(fac(w_ix))
                end do
            end if

            do w_ix=1, CTrans%NumSources - State%num_extra_redshiftwindows - num_cmb_freq*2
                do w_ix2=w_ix,CTrans%NumSources - State%num_extra_redshiftwindows - num_cmb_freq*2
                    iCl_Array(j,w_ix,w_ix2) =iCl_Array(j,w_ix,w_ix2) &
                        *fac(w_ix)*fac(w_ix2)*dbletmp
                    iCl_Array(j,w_ix2,w_ix) = iCl_Array(j,w_ix,w_ix2)
                end do
            end do
        end if

        iCl_scalar(j,C_Temp)  =  iCl_scalar(j,C_Temp)*dbletmp
        iCl_scalar(j,C_E) =  iCl_scalar(j,C_E)*dbletmp*ctnorm
        iCl_scalar(j,C_Cross) =  iCl_scalar(j,C_Cross)*dbletmp*sqrt(ctnorm)
        if (CTrans%NumSources>2) then
            iCl_scalar(j,C_Phi) = CP%ALens*iCl_scalar(j,C_Phi)*const_fourpi*ell**4
            !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
            !We put pix extra factors of l here to improve interpolation in l
            iCl_scalar(j,C_PhiTemp) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiTemp)*const_fourpi*ell**3
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
            iCl_scalar(j,C_PhiE) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiE)*const_fourpi*ell**3*sqrt(ctnorm)
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi E}
        end if
    end do
#ifndef __INTEL_COMPILER
    !$OMP END PARALLEL DO
#endif

    end subroutine CalcScalCls
