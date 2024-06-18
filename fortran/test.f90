
            if (.not. DoInt .or. UseLimber(ThisCT%ls%l(j)) .and. CP%WantScalars) then
                !Limber approximation for small scale lensing (better than poor version of above integral)
                xf = State%tau0-(ThisCT%ls%l(j)+0.5_dl)/IV%q
                if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
                    n=State%TimeSteps%IndexOf(xf)
                    xf= (xf-State%TimeSteps%points(n))/(State%TimeSteps%points(n+1)-State%TimeSteps%points(n))
                    sums(3) = (IV%Source_q(n,3)*(1-xf) + xf*IV%Source_q(n+1,3))*&
                        sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/IV%q
                else
                    sums(3)=0
                end if
            end if
            if (.not. DoInt .and. ThisSources%NonCustomSourceNum>3) then
                if (any(ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum)==0 .or. &
                    ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum) > j)) then
                    !When CMB does not need integral but other sources do
                    do n= State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows), &
                        min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))
                        !Full Bessel integration
                        a2 = aa(n)
                        bes_ix = bes_index(n)

                        J_l = a2 * ajl(bes_ix, j) + (1 - a2) * (ajl(bes_ix + 1, j) -&
                            ((a2 + 1) * ajlpr(bes_ix, j) + (2 - a2) * &
                            ajlpr(bes_ix + 1, j)) * fac(n)) !cubic spline
                        J_l = J_l * State%TimeSteps%dpoints(n)

                        sums(4) = sums(4) + IV%Source_q(n, 4) * J_l
                        do s_ix = 5, ThisSources%NonCustomSourceNum
                            sums(s_ix) = sums(s_ix) + IV%Source_q(n, s_ix) * J_l
                        end do
                    end do
                end if
            end if
        end if
