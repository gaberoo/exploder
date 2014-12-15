      subroutine initmod (odeparms)
        external odeparms
        double precision parms(8)
        common /myparms/ parms
        call odeparms(8,parms)
        return
      end

      subroutine SIS (I, t)
        implicit none
        double precision, intent(in)  :: t
        double precision, intent(out) :: I
        double precision :: a, b
        double precision :: N, beta, mu, psi, rho, tmax, scond, I0
        common /myparms/ N, beta, mu, psi, rho, tmax, scond, I0
        if (N .gt. 0.0) then
          a = beta - (mu+psi)
          b = beta/N
          I = a/(b+(a-b)*exp(-a*t))
        else
          a = beta - (mu+psi)
          I = exp(a*t)
        end if
        return
      end

      subroutine logGrowth (I, t)
        implicit none
        double precision, intent(out) :: I
        double precision, intent(in) :: t
        double precision :: a, b, R0
        double precision :: N, beta, mu, psi, rho, tmax, scond, I0
        common /myparms/ N, beta, mu, psi, rho, tmax, scond, I0
        if (N .gt. 0.0) then
          a = beta - (mu+psi)
          b = beta/N
          I = a/(b+(a/I0-b)*exp(-a*t))
        else
          a = beta - (mu+psi)
          I = exp(a*t)
        end if
      end

      subroutine derivs (neq, t, y, ydot, yout, ip)
        implicit none
        double precision :: t
        double precision, dimension(3) :: y, ydot
        double precision, dimension(*) :: yout
        double precision :: N, beta, mu, psi, rho, tmax, scond, I0
        double precision :: I, fSI, dsamp, dcoal
        integer neq, ip(*)
        common /myparms/ N, beta, mu, psi, rho, tmax, scond, I0
        if (ip(1) < 1) call rexit("nout should be at least 1")
c       call SIS(I,tmax-t)
        call logGrowth(I,tmax-t)
        if (N .gt. 0.0) then
          fSI = beta*(1-I/N)*I
        else
          fSI = beta*I
        end if
        dcoal = -fSI*(y(1)/I)**2
        dsamp = psi*I
        if (scond .gt. 0.0) then
          ydot(1) = dcoal
        else
          ydot(1) = dsamp + dcoal
        end if
        ydot(2) = dsamp
        ydot(3) = dcoal
        yout(1) = dsamp
        yout(2) = dcoal
        yout(3) = I
        return
      end

      subroutine derivsLI (neq, t, y, ydot, yout, ip)
        implicit none
        double precision :: t
        double precision, dimension(4) :: y, ydot
        double precision, dimension(*) :: yout
        double precision :: N, beta, mu, psi, rho, tmax, scond, I0
        double precision :: I, fSI, dsamp, dcoal
        integer :: neq, ip(*)
        common /myparms/ N, beta, mu, psi, rho, tmax, scond, I0
        if (ip(1) < 1) call rexit("nout should be at least 1")
        I = y(4)
        if (N .gt. 0.0) then
          fSI = beta*(1-y(4)/N)*y(4)
        else
          fSI = beta*y(4)
        end if
        dcoal = -fSI*(y(1)/y(4))**2
        dsamp = psi*y(4)
        if (scond .gt. 0.0) then
          ydot(1) = dcoal
        else
          ydot(1) = dsamp + dcoal
        end if
        ydot(2) = dsamp
        ydot(3) = dcoal
        ydot(4) = -fSI + mu*y(4)
        yout(1) = dsamp
        yout(2) = dcoal
        return
      end

      subroutine jac (neq, t, y, ml, mmu, pd, nrowpd, yout, ip)
        integer :: neq, ml, mmu, nrowpd, ip
        double precision :: y(*), pd(nrowpd,*), yout(*), t
        double precision :: N, beta, mu, psi, rho, tmax, scond, I0
        double precision :: I, fSI
        common /myparms/ N, beta, mu, psi, rho, tmax, scond, I0
        call SIS(I,tmax-t)
        if (N .gt. 0.0) then
          fSI = beta*(1-I/N)*I
        else
          fSI = beta*I
        end if
        if (N .gt. 0.0) then
c         if (y(1) .gt. 1.0 .and. I .gt. 1.0) then
c           pd(1,1) = -psi-fSI/(I*(I-1))*(2*y(1)-1)
            pd(1,1) = -psi-2*fSI*y(1)/I**2
c         else
c           pd(1,1) = 0.0
c         end if
        else
c         if (y(1) .gt. 1.0 .and. I .gt. 1.0) then
c           pd(1,1) = -psi-fSI/I**2*(2*y(1)-1)
            pd(1,1) = -psi-2*fSI*y(1)/I**2
c         else
c           pd(1,1) = 0.0
c         end if
        end if
        pd(1,2) = 0.0
        pd(1,3) = 0.0
        pd(2,1) = -psi
        pd(2,2) = 0.0
        pd(2,3) = 0.0
        pd(3,1) = pd(1,1)+psi
        pd(3,2) = 0.0
        pd(3,3) = 0.0
        return 
      end


