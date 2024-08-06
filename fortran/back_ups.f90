subroutine derivs(EV,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the scalar perturbations
    use constants, only : barssc0, Compton_CT, line21_const
    use MassiveNu
    use Recombination
    implicit none
    type(EvolutionVars) EV
    integer n,nu_i
    ! andrea
    real(dl) ay(n),ayprime(n), ayin(n), ayprimein(n)
    real(dl) drhophot,piphot, qgphot,opac_qgphot, opac_phot, octgprime
    real(dl) polter_freq, opac_rayleigh, opac_tot
    integer f_i, ix_off, s_ix
    ! andrea
    real(dl) tau, w
    real(dl) k,k2
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t,sigma,polter
    real(dl) w_dark_energy_t !equation of state of dark energy
    real(dl) gpres_noDE !Pressure with matter and radiation, no dark energy
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) E2
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz
    real(dl) dgpi,dgrho_matter,grho_matter, clxnu, gpres_nu
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    real(dl) xe,Trad, Delta_TM, Tmat, Delta_TCMB
    real(dl) delta_p_b, wing_t, wing2_t,winv_t
    real(dl) Delta_source2, polter_line
    real(dl) Delta_xe, Tspin, tau_eps, tau_fac, Tb
    integer lineoff,lineoffpol
    !Variables for source calculation
    real(dl) diff_rhopi, pidot_sum, dgpi_diff, phi
    real(dl) E(2:3), Edot(2:3)
    real(dl) phidot, polterdot, polterddot, octg, octgdot
    real(dl) lenswindow
    real(dl) :: opacity(nscatter), dopacity(nscatter), &
                            ddopacity(nscatter), exptau(nscatter), &
                            visibility(nscatter), dvisibility(nscatter), &
                            ddvisibility(nscatter)
    real(dl) ISW, quadrupole_source, doppler, monopole_source, tau0, ang_dist
    real(dl) dgrho_de, dgq_de, cs2_de

    k=EV%k_buf
    k2=EV%k2_buf

    !  Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call EV%ThermoData%Values_array(tau,a,cs2,opacity,dopacity)
    else
        call EV%ThermoData%Values_array(tau,a,cs2,opacity)
    end if
    a2=a*a

    etak=ay(ix_etak)

    !  CDM variables
    clxc=ay(ix_clxc)

    !  Baryon variables
    clxb=ay(ix_clxb)
    vb=ay(ix_vb)
    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2

    if (EV%is_cosmological_constant) then
        grhov_t = State%grhov * a2
        w_dark_energy_t = -1_dl
    else
        call State%CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)
    end if

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    gpres_nu=0
    grhonu_t=0

    if (State%CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,grhonu_t,gpres_nu,dgrho_matter,dgq, wnu_arr)
    end if

    grho_matter=grhonu_t+grhob_t+grhoc_t
    grho = grho_matter+grhor_t+grhog_t+grhov_t
    gpres_noDE = gpres_nu + (grhor_t + grhog_t)/3

    if (State%flat) then
        adotoa=sqrt(grho/3)
        cothxor=1._dl/tau
    else
        adotoa=sqrt((grho+State%grhok)/3._dl)
        cothxor=1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
    end if

    dgrho = dgrho_matter

    if (EV%no_nu_multpoles) then
        !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
        !Approximate total density variables with just matter terms
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxr=-4*dz/k
        qr=-4._dl/3*z
        pir=0
    else
        !  Massless neutrinos
        clxr=ay(EV%r_ix)
        qr  =ay(EV%r_ix+1)
        pir =ay(EV%r_ix+2)
    endif

    pig=0
    if (EV%no_phot_multpoles) then
        if (.not. EV%no_nu_multpoles) then
            z=(0.5_dl*dgrho/k + etak)/adotoa
            dz= -adotoa*z - 0.5_dl*dgrho/k
            clxg=-4*dz/k-4/k*opacity(1)*(vb+z)
            qg=-4._dl/3*z
        else
            clxg=clxr-4/k*opacity(1)*(vb+z)
            qg=qr
        end if
    else
        !  Photons
        clxg=ay(EV%g_ix)
        qg=ay(EV%g_ix+1)
        if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if

    !  8*pi*a*a*SUM[rho_i*clx_i] - radiation terms
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=dgq + grhog_t*qg+grhor_t*qr

    if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles )  then
        call Phot_Integrate_L012(EV, ay, drhophot, qgphot, piphot, opacity, opac_qgphot, opac_phot)
        if (rayleigh_back) then
            dgrho=dgrho + grhog_t*drhophot
            dgq=dgq + grhog_t*qgphot
        else
            opac_qgphot = 0
            opac_phot = opacity(1)
        end if
    end if

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    if (.not. EV%is_cosmological_constant) then
        call State%CP%DarkEnergy%PerturbedStressEnergy(dgrho_de, dgq_de, &
            a, dgq, dgrho, grho, grhov_t, w_dark_energy_t, gpres_noDE, etak, &
            adotoa, k, EV%Kf(1), ay, ayprime, EV%w_ix)
        dgrho = dgrho + dgrho_de
        dgq = dgq + dgq_de
    end if

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    if (State%flat) then
        !eta*k equation
        sigma=(z+1.5_dl*dgq/k2)
        ayprime(ix_etak)=0.5_dl*dgq
    else
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
        ayprime(ix_etak)=0.5_dl*dgq + State%curv*z
    end if

    if (.not. EV%is_cosmological_constant) &
        call State%CP%DarkEnergy%PerturbationEvolve(ayprime, w_dark_energy_t, &
        EV%w_ix, a, adotoa, k, z, ay)

    !  CDM equation of motion
    clxcdot=-k*z
    ayprime(ix_clxc)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(ix_clxb)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    !Sources
    if (EV%Evolve_baryon_cs) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            Tmat = State%CP%Recomb%T_m(a)
        else
            Tmat = State%CP%TCMB/a
        end if
        if (EV%Evolve_TM) then
            Delta_TM = ay(EV%Tg_ix)
        else
            Delta_TM = clxg/4
        end if
        delta_p_b = barssc0*(1._dl-0.75d0*State%CP%yhe+(1._dl-State%CP%yhe)*opacity(1)*a2/State%akthom)*Tmat*(clxb + delta_tm)
    else
        Delta_TM = clxg/4
        delta_p_b = cs2*clxb
    end if

    if (State%CP%Evolve_delta_xe) then
        if (EV%saha) then
            xe=State%CP%Recomb%x_e(a)
            Delta_xe = (1-xe)/(2-xe)*(-clxb + (3._dl/2+  CB1/Tmat)*Delta_TM)
        else
            Delta_xe = ay(EV%xe_ix)
        end if
    else
        Delta_xe = 0
    end if

    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
        !  ddota/a
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        adotdota=(adotoa*adotoa-gpres)/2

        pig = 32._dl/45/opacity(1)*k*(sigma+vb)

        !  First-order approximation to baryon-photon splip
       slip = - (2*adotoa/(1+pb43) + dopacity(1)/opacity(1))* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity(1)*(1+pb43))

        if (second_order_tightcoupling) then
            ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming flat
            !AL: First order slip seems to be fine here to 2e-4

            !  8*pi*G*a*a*SUM[rho_i*sigma_i]
            dgs = grhog_t*pig+grhor_t*pir

            ! Define shear derivative to first order
            sigmadot = -2*adotoa*sigma-dgs/k+etak

            !Once know slip, recompute qgdot, pig, pigdot
            qgdot = k*(clxg/4._dl-pig/2._dl) +opacity(1)*slip

            pig = 32._dl/45/opacity(1)*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity(1)*11._dl/6._dl/opacity(1)**2)) &
                + (32._dl/45._dl/opacity(1)**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

            pigdot = -(32._dl/45._dl)*(dopacity(1)/opacity(1)**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
                dopacity(1)*11._dl/6._dl/opacity(1)**2 ) &
                + (32._dl/45._dl/opacity(1))*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
                *(dopacity(1)/opacity(1)**2))

            EV%pigdot = pigdot

        end if

        !  Use tight-coupling approximation for vb
        !  zeroth order approximation to vbdot + the pig term
        vbdot=(-adotoa*vb+cs2*k*clxb + k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip
        EV%pig = pig
    ! and
    else
        if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles)  then 
            vbdot = -adotoa*vb+cs2*k*clxb -photbar*(opac_phot*(4._dl/3*vb-qg) -opac_qgphot)
        else
            vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity(1)*(4._dl/3*vb-qg)
        end if
    end if
    ! and
    ayprime(ix_vb)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        !qgdot=4._dl/3*(-vbdot-adotoa*vb+delta_p_b*k)/pb43 &
        if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles)  then 
        qgdot=4._dl/3*(photbar*opacity(1)*(4._dl/3*vb-qg))/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
        else
        qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig  
        end if
        ayprime(EV%g_ix+1)=qgdot

        !  Use explicit equations for photon moments if appropriate
        if (.not. EV%tightcoupling) then
            E2=ay(EV%polind+2)
            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
            ix= EV%g_ix+2
            if (EV%lmaxg>2) then
                pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity(1)*(pig - polter) &
                    +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
                do  l=3,EV%lmaxg-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity(1)*ay(ix)
                end do
                ix=ix+1
                !  Truncate the photon moment expansion
                ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity(1)*ay(ix)
            else !closed case
                pigdot=EV%denlk(2)*qg-opacity(1)*(pig - polter) +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
            endif
            !  Polarization
            !l=2
            ix=EV%polind+2
            if (EV%lmaxgpol>2) then
                ayprime(ix) = -opacity(1)*(ay(ix) - polter) - k/3._dl*ay(ix+1)
                do l=3,EV%lmaxgpol-1
                    ix=ix+1
                    ayprime(ix)=-opacity(1)*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                end do
                ix=ix+1
                !truncate
                ayprime(ix)=-opacity(1)*ay(ix) + &
                    k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
            else !closed case
                ayprime(ix) = -opacity(1)*(ay(ix) - polter)
            endif
        end if
    end if
    !write(*,*) 'ray condi', EV%Rayleigh
    if (EV%Rayleigh) then
    !     !assume after tight coupling ended
         do f_i = 1, num_cmb_freq
         !       !write(*,*) 'freq_factors = ', freq_factors
                opac_rayleigh = State%CP%Recomb%Recombination_rayleigh_eff(a)*State%akthom/a2*(min(1._dl,&
                                freq_factors(f_i,1)/a2**2 + freq_factors(f_i,2)/a2**3  + freq_factors(f_i,3)/a2**4))
                opac_tot=opac_rayleigh+opacity(1)
                ind = EV%g_ix_freq + (f_i-1)*EV%freq_neq
                ayprime(ind)=-k*ay(ind+1)
                ayprime(ind+1)= EV%denlk(1)*ay(ind)-EV%denlk2(1)*ay(ind+2) + &
                                4._dl/3*photbar*(-opacity(1)*ay(ind+1) + opac_rayleigh*(4._dl/3*vb-qg-ay(ind+1)))/pb43
                E2=ay(EV%polind_freq + (f_i-1)*EV%freq_neq +2)
                polter_freq = ay(ind+2)/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
                if (EV%lmaxg>2) then
                     ayprime(ind+2)= EV%denlk(2)*ay(ind+1)-EV%denlk2(2)*ay(ind+3)-opac_tot*(ay(ind+2) - polter_freq) &
                                    -opac_rayleigh*(ay(EV%g_ix+2) - polter) 
                     ix = ind+2
                     do  l=3,EV%lmaxg-1
                        ix = ix+1
                        ayprime(ix) = (EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opac_tot*ay(ix)-(opac_rayleigh)*ay(EV%g_ix+l)
                     end do
                     ix=ix+1
                  !  Truncate the photon moment expansion
                     ayprime(ix) = k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opac_tot*ay(ix)-opac_rayleigh*ay(EV%g_ix+EV%lmaxg)
                else 
                     ayprime(ind+2) = EV%denlk(2)*ay(ind+1)-opac_tot*(ay(ind+2) - polter_freq)-opac_rayleigh*(ay(EV%g_ix+2) - polter) 
                endif
        !  Polarization
                l=2
                ix=EV%polind_freq+(f_i-1)*EV%freq_neq + 2
                if (EV%lmaxgpol>2) then
                    ayprime(ix) = -opac_tot*(ay(ix) - polter_freq) - k/3._dl*ay(ix+1)- (opac_rayleigh)*(ay(EV%polind+2) - polter)
                    do l=3,EV%lmaxgpol-1
                        ix = ix+1
                        ayprime(ix) = -opac_tot*ay(ix)-opac_rayleigh*ay(EV%polind+l) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                    end do
                    ix=ix+1
                    !truncate
                    ayprime(ix) = -opac_tot*ay(ix)-opac_rayleigh*ay(EV%polind+EV%lmaxgpol) + &
                                k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
                else !closed case
                    ayprime(ix) = -opac_tot*(ay(ix) - polter_freq)-opac_rayleigh*ay(EV%polind+2)
                endif
         end do
     end if

    if (.not. EV%no_nu_multpoles) then
        !  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%r_ix)=clxrdot
        qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
        ayprime(EV%r_ix+1)=qrdot
        if (EV%high_ktau_neutrino_approx) then
            !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
            !Method from arXiv:1104.2933
            !                if (.not. EV%TightCoupling) then
            !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            !                 adotdota=(adotoa*adotoa-gpres)/2
            !                end if
            !                ddz=(2*adotoa**2 - adotdota)*z  &
            !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*(grhoc_t*clxc+grhob_t*clxb) ) &
            !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + grhoc_t*clxcdot + grhob_t*clxbdot )
            !                dz= -adotoa*z - 0.5_dl*dgrho/k
            !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
            pirdot= -3*pir*cothxor - clxrdot
            ayprime(EV%r_ix+2)=pirdot

            !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
            !                ayprime(EV%lmaxg+9)=pirdot
            !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
            !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
            !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
        else
            ix=EV%r_ix+2
            if (EV%lmaxnr>2) then
                pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
                do l=3,EV%lmaxnr-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
                end do
                !  Truncate the neutrino expansion
                ix=ix+1
                ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
            else
                pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
            end if
        end if
    end if ! no_nu_multpoles

    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            Delta_TCMB = clxg/4
            xe = State%CP%Recomb%x_e(a)
            Trad = State%CP%TCMB/a

            !Matter temperature
            !Recfast_CT = (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc [a_R = radiation constant]
            ayprime(EV%Tg_ix) = -2*k*(z+vb)/3 - a*  Compton_CT * (Trad**4) * xe / (1._dl+xe+State%fHe) * &
                ((1- Trad/Tmat)*(Delta_TCMB*4 + Delta_xe/(1+xe/(1+State%fHe))) &
                + Trad/Tmat*(Delta_Tm - Delta_TCMB)  )

            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) =  Get21cm_dTs(a,clxb,ay(EV%Ts_ix),Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
            end if
        else
            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) = -k*(4._dl/3._dl*z+qg)/4  !Assume follows Delta_TM which follows clxg
            end if
        end if
    end if

    if (State%CP%Evolve_delta_xe .and. .not. EV%saha) then
        ayprime(EV%xe_ix) = &
            State%CP%Recomb%dDeltaxe_dtau(a, Delta_xe,clxb, Delta_Tm, k*z/3,k*vb, adotoa)
    end if

    if (State%CP%Do21cm) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            if (State%CP%SourceTerms%line_reionization) then
                lineoff = EV%reion_line_ix+1
                lineoffpol = lineoff+EV%lmaxline-1

                if (tau> EV%ThermoData%tau_start_redshiftwindows) then
                    !Multipoles of 21cm

                    polter_line = ay(lineoff+2)/10+9._dl/15*ay(lineoffpol+2)

                    call interp_window(State%TimeSteps,State%Redshift_W(1),tau,wing_t,wing2_t,winv_t)

                    delta_source2 = Get21cm_source2(a,clxb,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe,k*(z+vb)/adotoa/3)


                    !Drop some small terms since mulipoles only enter into reionzation anyway
                    !monopole
                    ayprime(lineoff) = -k*ay(lineoff+1) +  wing_t * clxb + wing2_t*delta_source2 + k*z/3*winV_t


                    !dipole
                    ayprime(lineoff+1)= EV%denlk(1)*ay(lineoff)-EV%denlk2(1)*ay(lineoff+2) - opacity(1)*ay(lineoff+1) &
                        -wing2_t * ( qg/4 - vb/3)   ! vb/3*WinV_t)

                    !quadrupole
                    ayprime(lineoff+2)= EV%denlk(2)*ay(lineoff+1)-EV%denlk2(2)*ay(lineoff+3) &
                        +opacity(1)*(polter_line -ay(lineoff+2) ) -   2._dl/15*k*sigma*winV_t &
                        - wing2_t * ay(EV%g_ix+2)/4

                    do  l=3,EV%lmaxline-1
                        ayprime(lineoff+l)=EV%denlk(l)*ay(lineoff+l-1)-EV%denlk2(l)*ay(lineoff+l+1)-opacity(1)*ay(lineoff+l) &
                            - wing2_t * ay(EV%g_ix+l)/4
                    end do
                    !truncate
                    ayprime(lineoff+EV%lmaxline)=k*ay(lineoff+EV%lmaxline-1)-(EV%lmaxline+1)*cothxor*ay(lineoff+EV%lmaxline)  &
                        -opacity(1)*ay(lineoff+EV%lmaxline) - wing2_t * ay(EV%g_ix+EV%lmaxline)/4

                    !  21cm Polarization
                    !l=2
                    ayprime(lineoffpol+2) = -opacity(1)*(ay(lineoffpol+2) - polter_line) - k/3._dl*ay(lineoffpol+3)
                    !and the rest
                    do l=3,EV%lmaxline-1
                        ayprime(lineoffpol+l)=-opacity(1)*ay(lineoffpol+l) + EV%denlk(l)*ay(lineoffpol+l-1) -&
                            EV%polfack(l)*ay(lineoffpol+l+1)
                    end do

                    !truncate
                    ayprime(lineoffpol+EV%lmaxline)=-opacity(1)*ay(lineoffpol+EV%lmaxline) + &
                        k*EV%poltruncfac*ay(lineoffpol+EV%lmaxline-1)-(EV%lmaxline+3)*cothxor*ay(lineoffpol+EV%lmaxline)
                else
                    ayprime(lineoff:lineoffpol+EV%lmaxline)=0
                end if
            end if
        end if
    end if



    !  Massive neutrino equations of motion.
    if (State%CP%Num_Nu_massive >0) then
        !DIR$ LOOP COUNT MIN(1), AVG(1)
        do nu_i = 1, State%CP%Nu_mass_eigenstates
            if (EV%MassiveNuApprox(nu_i)) then
                !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
                !see astro-ph/0203507
                G11_t=EV%G11(nu_i)/a/a2
                G30_t=EV%G30(nu_i)/a/a2
                off_ix = EV%nu_ix(nu_i)
                w=wnu_arr(nu_i)
                ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
                ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
                ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
                ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
            else
                ind=EV%nu_ix(nu_i)
                !DIR$ LOOP COUNT MIN(3), AVG(3)
                do i=1,EV%nq(nu_i)
                    q=State%NuPerturbations%nu_q(i)
                    aq=a*State%nu_masses(nu_i)/q
                    v=1._dl/sqrt(1._dl+aq*aq)

                    ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
                    ind=ind+1
                    ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                    ind=ind+1
                    if (EV%lmaxnu_tau(nu_i)==2) then
                        ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                    else
                        ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                            +k*8._dl/15._dl*sigma
                        do l=3,EV%lmaxnu_tau(nu_i)-1
                            ind=ind+1
                            ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                        end do
                        !  Truncate moment expansion.
                        ind = ind+1
                        ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                    end if
                    ind = ind+1
                end do
            end if
        end do

        if (EV%has_nu_relativistic) then
            ind=EV%nu_pert_ix
            ayprime(ind)=+k*a2*qr -k*ay(ind+1)
            ind2= EV%r_ix
            do l=1,EV%lmaxnu_pert-1
                ind=ind+1
                ind2=ind2+1
                ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                    +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
            end do
            ind=ind+1
            ind2=ind2+1
            ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
        end if
    end if
    ! and
    ayprimein = ayprime
    ! and
    if (associated(EV%OutputTransfer) .or. associated(EV%OutputSources)) then
        if (EV%TightCoupling .or. EV%no_phot_multpoles) then
            E=0
            Edot=0
        else
            E = ay(EV%polind+2:EV%polind+3)
            Edot = ayprime(EV%polind+2:EV%polind+3)
        end if
        if (EV%no_nu_multpoles) then
            pirdot=0
            qrdot = -4*dz/3
        end if
        if (EV%no_phot_multpoles) then
            pigdot=0
            octg=0
            octgdot=0
            qgdot = -4*dz/3
        else
            if (EV%TightCoupling) then
                if (second_order_tightcoupling) then
                    octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity(1))
                    E(2) = pig/4 + pigdot*(1._dl/opacity(1))*(-5._dl/8._dl)
                    E(3) = (3._dl/7._dl)*(EV%k_buf/opacity(1))*E(2)
                    Edot(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity(1)/opacity(1)**2))
                else
                    pigdot = -dopacity(1)/opacity(1)*pig + 32._dl/45*k/opacity(1)*(-2*adotoa*sigma  &
                        +etak/EV%Kf(1)-  dgpi/k +vbdot )
                    Edot(2) = pigdot/4
                    E(2) = pig/4
                    octg=0
                end if
                octgdot=0
            else
                octg=ay(EV%g_ix+3)
                octgdot=ayprime(EV%g_ix+3)
            end if
        end if
        if (EV%is_cosmological_constant) then
            dgrho_de=0
            dgq_de=0
        end if

        dgpi  = grhor_t*pir + grhog_t*pig
        dgpi_diff = 0  !sum (3*p_nu -rho_nu)*pi_nu
        pidot_sum = grhog_t*pigdot + grhor_t*pirdot
        clxnu =0
        if (State%CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, adotoa, dgpi=dgpi, clxnu_all=clxnu, &
                dgpi_diff=dgpi_diff, pidot_sum=pidot_sum)
        end if
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff)*adotoa + &
            State%CP%DarkEnergy%diff_rhopi_Add_Term(dgrho_de, dgq_de, grho, &
            gpres, w_dark_energy_t, State%grhok, adotoa, &
            EV%kf(1), k, grhov_t, z, k2, ayprime, ay, EV%w_ix)
        phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)

        if (associated(EV%OutputTransfer)) then
            EV%OutputTransfer(Transfer_kh) = k/(State%CP%h0/100._dl)
            EV%OutputTransfer(Transfer_cdm) = clxc
            EV%OutputTransfer(Transfer_b) = clxb
            EV%OutputTransfer(Transfer_g) = clxg
            EV%OutputTransfer(Transfer_r) = clxr
            EV%OutputTransfer(Transfer_nu) = clxnu
            EV%OutputTransfer(Transfer_tot) =  dgrho_matter/grho_matter !includes neutrinos
            EV%OutputTransfer(Transfer_nonu) = (grhob_t*clxb+grhoc_t*clxc)/(grhob_t + grhoc_t)
            EV%OutputTransfer(Transfer_tot_de) =  dgrho/grho_matter
            !Transfer_Weyl is k^2Phi, where Phi is the Weyl potential
            EV%OutputTransfer(Transfer_Weyl) = k2*phi
            EV%OutputTransfer(Transfer_Newt_vel_cdm)=  -k*sigma/adotoa
            EV%OutputTransfer(Transfer_Newt_vel_baryon) = -k*(vb + sigma)/adotoa
            EV%OutputTransfer(Transfer_vel_baryon_cdm) = vb
            if (State%CP%do21cm) then
                Tspin = State%CP%Recomb%T_s(a)
                xe = State%CP%Recomb%x_e(a)

                tau_eps = a*line21_const*State%NNow/a**3/adotoa/Tspin/1000
                delta_source2 = Get21cm_source2(a,clxb,clxg/4,Delta_Tm,Delta_xe,Tmat,&
                    State%CP%TCMB/a,xe,k*(z+vb)/adotoa/3)
                tau_fac = tau_eps/(exp(tau_eps)-1)
                EV%OutputTransfer(Transfer_monopole) = ( clxb + Trad/(Tspin-Trad)*delta_source2 )  &
                    + (tau_fac-1)*(clxb - (delta_source2 + clxg/4)  )

                EV%OutputTransfer(Transfer_vnewt) = tau_fac*k*(vb+sigma)/adotoa
                EV%OutputTransfer(Transfer_Tmat) =  delta_TM
                if (State%CP%SourceTerms%use_21cm_mK) then
                    Tb = (1-exp(-tau_eps))*a*(Tspin-Trad)*1000

                    EV%OutputTransfer(Transfer_monopole) = EV%OutputTransfer(Transfer_monopole)*Tb
                    EV%OutputTransfer(Transfer_vnewt) = EV%OutputTransfer(Transfer_vnewt)*Tb
                    EV%OutputTransfer(Transfer_Tmat) = EV%OutputTransfer(Transfer_Tmat)*Tb
                end if
            end if
        end if
        if (associated(EV%OutputSources)) then

            EV%OutputSources = 0
            call EV%ThermoData%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
                visibility, dvisibility, ddvisibility, exptau, lenswindow)

            tau0 = State%tau0
            phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                diff_rhopi+ k*sigma*(gpres + grho))/k2
            !time derivative of shear
            sigmadot = -adotoa*sigma - 1.0d0/2.0d0*dgpi/k + k*phi
            !quadrupole source derivatives; polter = pi_g/10 + 3/5 E_2
            polter = pig/10+9._dl/15*E(2)
            polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                         k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                         50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                         *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                         polter - 3.0d0/5.0d0*E(2))*dopacity(1) + (-1.0d0/10.0d0*pigdot &
                        + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity(1)
            !Temperature source terms, after integrating by parts in conformal time

            s_ix = 0
            do f_i = 1, nscatter
                if (.not. EV%no_phot_multpoles .and. EV%Rayleigh .and. f_i>1) then
                    ix_off=EV%g_ix_freq+(f_i-2)*EV%freq_neq
                    ay(EV%g_ix:EV%g_ix+EV%freq_neq-1) = ayin(EV%g_ix:EV%g_ix+EV%freq_neq-1)+ &
                                                        ayin(ix_off:ix_off+EV%freq_neq-1)
                    ayprime(EV%g_ix:EV%g_ix+EV%freq_neq-1) = ayprimein(EV%g_ix:EV%g_ix+EV%freq_neq-1)+ &
                    ayprimein(ix_off:ix_off+EV%freq_neq-1)
                    pig = ay(EV%g_ix+2)
                    pigdot = ayprime(EV%g_ix+2)
                    octg = ay(EV%g_ix+3)
                    octgprime = ayprime(EV%g_ix+3)
                    clxg = ay(EV%g_ix)
                    qg  = ay(EV%g_ix+1)
                    qgdot = ayprime(EV%g_ix+1)
                    polter = pig/10+9._dl/15*E(2)
                    polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
                    polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                         k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                         50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                         *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                         polter - 3.0d0/5.0d0*E(2))*dopacity(f_i) + (-1.0d0/10.0d0*pigdot &
                        + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity(f_i)
                end if
                !2phi' term (\phi' + \psi' in Newtonian gauge), phi is the Weyl potential
                ISW = 2*phidot*exptau(f_i)
                monopole_source =  (-etak/(k*EV%Kf(1)) + 2*phi + clxg/4)*visibility(f_i)
                doppler = ((sigma + vb)*dvisibility(f_i) + (sigmadot + vbdot)*visibility(f_i))/k
                quadrupole_source = (5.0d0/8.0d0)*(3*polter*ddvisibility(f_i) + 6*polterdot*dvisibility(f_i) &
                                    + (k**2*polter + 3*polterddot)*visibility(f_i))/k**2

                EV%OutputSources(s_ix+1) = ISW + doppler + monopole_source + quadrupole_source
                ang_dist = f_K(tau0-tau)

                if (tau < tau0) then
                    !E polarization source
                    EV%OutputSources(s_ix+2)= visibility(f_i)*polter*(15._dl/8._dl)/(ang_dist**2*k2)
                    !factor of four because no 1/16 later
                end if

                if (rayleigh_diff .and. f_i>2) then
                    EV%OutputSources(s_ix+1:s_ix+2)= EV%OutputSources(s_ix+1:s_ix+2) - EV%OutputSources(1:2)
                end if
                s_ix = s_ix+2
            
                if (size(EV%OutputSources) > 2 .and. f_i==1) then
                    ! and
                    !Get lensing sources
                    if (tau>State%tau_maxvis .and. tau0-tau > 0.1_dl) then
                        EV%OutputSources(3) = -2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
                        !write(*,*) 's3 : ', 2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
                        !We include the lensing factor of two here
                    end if
                    s_ix=s_ix+1
                end if
                !EV%OutputSources(s_ix) = 0
            end do !f_1
        end if
    end if
    
    end subroutine derivs
