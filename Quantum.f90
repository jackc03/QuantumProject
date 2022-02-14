Program Quantum
    implicit none
    complex*16                  :: rho(2,2)
    double precision            :: x(2,2)
    double precision            :: t, tf, ti
    double precision            :: tw = 2d-14
    double precision            :: w
    double precision            :: dt, E2
    double precision, parameter :: eV = 1.602d-19
    double precision, parameter :: hbar = 1.05d-34
    double precision, parameter :: pi =3.1415926535897932
    integer                     :: i, Nt
    complex*16                  :: electricfield(2,2)

    rho(1,1) = 1
    rho(1,2) = 0
    rho(2,1) = 0
    rho(2,2) = 0
    E2 = 1.5 * eV
    w = E2/hbar 
    dt = (2*pi)/(100 * w)
    ti = -100d-15
    tf = 100d-15
    Nt = 7330
    t = ti
    x(1,2) = 10d-10
    x(2,1) = x(1,2)

    call openfiles()

    do i = 1, Nt
        call RK4(rho, t, dt)
        call Printfiles(rho,t)
        write(105,*) t, Et(t)
        write(106,*) t, real(Trace(rho,-eV*x)) * (1d28) ! 10^28 atoms per m^3
        t = t + dt
    end do

    call closefiles()
    
contains

    subroutine RK4(rho, t, dt)
        complex*16, intent(inout)       :: rho(2,2)
        double precision, intent(inout) :: t
        double precision, intent(in)    :: dt
        complex*16                      :: K1(2,2), K2(2,2), K3(2,2), K4(2,2) 
        
        K1 = rhoDOT(rho         , t         ) * dt
        K2 = rhoDOT(rho + K1/2d0, t + dt/2d0) * dt
        K3 = rhoDOT(rho + K2/2d0, t + dt/2d0) * dt
        K4 = rhoDOT(rho + K3    , t + dt    ) * dt


        rho = rho + (1/6d0) * (K1 + 2 * (K2+K3) + K4)
    end subroutine

    
    function rhoDOT(rho, t)
        complex*16                  :: rhoDOT(2,2)
        complex*16, intent(in)      :: rho(2,2)
        double precision            :: x(2,2)
        double precision,intent(in) :: t
        double precision            :: E1,E2 
        complex*16                  :: ii = (0d0, 1d0)
        double precision            :: e0 = 1.602d-19
        x(1,2) = 10d-10
        x(2,1) = x(1,2)
        E1   = 0
        E2   = 1.5 * e0

        rhoDOT(1,1) = (1/(ii*hbar)) * (rho(2,1) * x(1,2) - x(2,1) * rho(1,2)) * ABS(e0) * Et(t)
        rhoDOT(1,2) = (1/(ii*hbar)) * (rho(1,2)*(E1-E2) + (ABS(e0) * Et(t) * x(1,2)) * (rho(2,2) - rho(1,1)))
        rhoDOT(2,1) = (1/(ii*hbar)) * (rho(2,1)*(E2-E1) + (ABS(e0) * Et(t) * x(2,1)) * (rho(1,1) - rho(2,2)))
        rhoDOT(2,2) = (1/(ii*hbar)) * (rho(1,2) * x(2,1) - x(1,2) * rho(2,1)) * ABS(e0) * Et(t)
    end function

    function Et(t)
        double precision :: t, w, Emax
        double precision :: Et
        character(100)   :: numchar
        double precision :: E2 = 1.5 * 1.602d-19
        double precision :: hbar = 1.05d-34
        double precision :: tw = 2d-14


        w = (E2/hbar)*1.0d0
        Emax = 10**8
        Et = Emax * cos(w * t)*exp(-(t**2)/(tw**2))
    end function

    function Trace(M1, M2)
        complex*16       :: M1(2,2), MultMatrix(2,2), Trace
        double precision :: M2(2,2)

        MultMatrix = matmul(M1, M2)
        Trace = MultMatrix(1,1) + MultMatrix(2,2)
    end function

    subroutine Printfiles(rho,t)
        double precision       :: t
        complex*16, intent(in) :: rho(2,2)
        write(101,*) t, real(rho(1,1)), aimag(rho(1,1))
        write(102,*) t, real(rho(1,2)), aimag(rho(1,2))
        write(103,*) t, real(rho(2,1)), aimag(rho(2,1))
        write(104,*) t, real(rho(2,2)), aimag(rho(2,2))
    end subroutine

    subroutine openfiles()
        open(file='RHO_11.dat', unit=101)
        open(file='RHO_12.dat', unit=102)
        open(file='RHO_21.dat', unit=103)
        open(file='RHO_22.dat', unit=104)
        open(file='E.dat', unit=105)
        open(file='P.dat', unit=106)
    end subroutine

    subroutine closefiles()
        close(101)
        close(102)
        close(103)
        close(104)
        close(105)
        close(106)
    end subroutine
end Program    