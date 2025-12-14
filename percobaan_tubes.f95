! ==============================================================================
! TUGAS BESAR OS2103: KOMPUTASI OSEANOGRAFI
! TOPIK: ALGORITMA UNESCO NO. 44 (1983)
! FULL IMPLEMENTATION (ALGO 1-9)
! ==============================================================================

MODULE UNESCO_MATH
    IMPLICIT NONE
    
    ! Konstanta Global (Paper Hal 12)
    REAL, PARAMETER :: k_sal = 0.0162

CONTAINS

    ! =========================================================================
    ! [cite_start]ALGO 1: CONDUCTIVITY TO SALINITY (PSS-78) [cite: 1638-1708]
    ! =========================================================================
    FUNCTION ALGO1_SAL78(R_in, T_in, P_in) RESULT(S_out)
        REAL, INTENT(IN) :: R_in, T_in, P_in
        REAL :: S_out, RT, Rp, rt_val, DT, R_t
        ! Konstanta (Paper Hal 12)
        REAL, PARAMETER :: a0=0.0080, a1=-0.1692, a2=25.3851, a3=14.0941, a4=-7.0261, a5=2.7081
        REAL, PARAMETER :: b0=0.0005, b1=-0.0056, b2=-0.0066, b3=-0.0375, b4=0.0636, b5=-0.0144
        REAL, PARAMETER :: c0=0.6766097, c1=2.00564E-2, c2=1.104259E-4, c3=-6.9698E-7, c4=1.0031E-9
        REAL, PARAMETER :: d1=3.426E-2, d2=4.464E-4, d3=4.215E-1, d4=-3.107E-3
        REAL, PARAMETER :: e1=2.070E-5, e2=-6.370E-10, e3=3.989E-15

        DT = T_in - 15.0
        Rp = 1.0 + (P_in * (e1 + e2*P_in + e3*P_in**2)) / &
             (1.0 + d1*T_in + d2*T_in**2 + (d3 + d4*T_in)*R_in)
        rt_val = c0 + T_in*(c1 + T_in*(c2 + T_in*(c3 + T_in*c4)))
        R_t = R_in / (Rp * rt_val)
        RT = SQRT(ABS(R_t)) 
        S_out = a0 + RT*(a1 + RT*(a2 + RT*(a3 + RT*(a4 + RT*a5)))) + &
                (DT / (1.0 + k_sal*DT)) * &
                (b0 + RT*(b1 + RT*(b2 + RT*(b3 + RT*(b4 + RT*b5)))))
        IF (S_out .LT. 0.0) S_out = 0.0
    END FUNCTION ALGO1_SAL78

    ! =========================================================================
    ! [cite_start]ALGO 2: SALINITY TO CONDUCTIVITY (Iterasi Newton-Raphson) [cite: 1709-1733]
    ! =========================================================================
    FUNCTION ALGO2_COND(S_in, T_in, P_in) RESULT(R_out)
        REAL, INTENT(IN) :: S_in, T_in, P_in
        REAL :: R_out, RT_guess, SI_calc, DELS, DERIV, DT
        REAL :: Rp_est, rt_val, A_pol, B_pol, C_pol, R_t_sq
        INTEGER :: i
        ! Konstanta PSS-78
        REAL, PARAMETER :: a1=-0.1692, a2=25.3851, a3=14.0941, a4=-7.0261, a5=2.7081
        REAL, PARAMETER :: b1=-0.0056, b2=-0.0066, b3=-0.0375, b4=0.0636, b5=-0.0144
        REAL, PARAMETER :: c0=0.6766097, c1=2.00564E-2, c2=1.104259E-4, c3=-6.9698E-7, c4=1.0031E-9
        REAL, PARAMETER :: d1=3.426E-2, d2=4.464E-4, d3=4.215E-1, d4=-3.107E-3
        REAL, PARAMETER :: e1=2.070E-5, e2=-6.370E-10, e3=3.989E-15

        DT = T_in - 15.0
        RT_guess = SQRT(S_in / 35.0) 
        DELS = 1.0
        i = 0

        DO WHILE (DELS .GT. 1.0E-5 .AND. i .LT. 20)
            i = i + 1
            SI_calc = 0.0080 + RT_guess*(a1 + RT_guess*(a2 + RT_guess*(a3 + RT_guess*(a4 + RT_guess*a5)))) + &
                      (DT/(1.0+k_sal*DT)) * &
                      (0.0005 + RT_guess*(b1 + RT_guess*(b2 + RT_guess*(b3 + RT_guess*(b4 + RT_guess*b5)))))
            DERIV = a1 + RT_guess*(2.*a2 + RT_guess*(3.*a3 + RT_guess*(4.*a4 + RT_guess*5.*a5))) + &
                    (DT/(1.0+k_sal*DT)) * &
                    (b1 + RT_guess*(2.*b2 + RT_guess*(3.*b3 + RT_guess*(4.*b4 + RT_guess*5.*b5))))
            RT_guess = RT_guess + (S_in - SI_calc) / DERIV
            DELS = ABS(S_in - SI_calc)
        END DO

        rt_val = c0 + T_in*(c1 + T_in*(c2 + T_in*(c3 + T_in*c4)))
        A_pol = d3 + d4*T_in
        B_pol = 1.0 + d1*T_in + d2*T_in**2
        C_pol = P_in * (e1 + e2*P_in + e3*P_in**2)
        
        R_t_sq = RT_guess**2
        R_out = R_t_sq * rt_val 
        Rp_est = 1.0 + C_pol / (B_pol + A_pol*R_out) 
        R_out = R_out * Rp_est
    END FUNCTION ALGO2_COND

    ! =========================================================================
    ! [cite_start]ALGO 3: DENSITY (EOS-80) [cite: 1866-1906]
    ! =========================================================================
    FUNCTION ALGO3_DENSITY(S, T, P) RESULT(Rho)
        REAL, INTENT(IN) :: S, T, P
        REAL :: Rho, P_bar, rho_w, K0, K_w, A, B, BulkMod
        P_bar = P / 10.0

        rho_w = 999.842594 + 6.793952E-2*T - 9.095290E-3*T**2 + &
                1.001685E-4*T**3 - 1.120083E-6*T**4 + 6.536332E-9*T**5

        K0 = rho_w + (0.824493 - 0.0040899*T + 7.6438E-5*T**2 - 8.2467E-7*T**3 + 5.3875E-9*T**4)*S + &
             (-0.00572466 + 1.0227E-4*T - 1.6546E-6*T**2)*S*SQRT(S) + 4.8314E-4*S**2

        K_w = 19652.21 + 148.4206*T - 2.327105*T**2 + 1.360477E-2*T**3 - 5.155288E-5*T**4
        BulkMod = K_w + (54.6746 - 0.603459*T + 1.09987E-2*T**2 - 6.1670E-5*T**3)*S + &
                  (0.07944 + 0.016483*T - 5.3009E-4*T**2)*S*SQRT(S)
        
        A = 3.239908 + 1.43713E-3*T + 1.16092E-4*T**2 - 5.77905E-7*T**3 + &
            (2.2838E-3 - 1.0981E-5*T - 1.6078E-6*T**2)*S + 1.91075E-4*S*SQRT(S)
        B = 8.50935E-5 - 6.12293E-6*T + 5.2787E-8*T**2 + (-9.9348E-7 + 2.0816E-8*T + 9.1697E-10*T**2)*S
        
        BulkMod = BulkMod + A*P_bar + B*P_bar**2
        Rho = K0 / (1.0 - P_bar / BulkMod)
    END FUNCTION ALGO3_DENSITY

    ! =========================================================================
    ! [cite_start]ALGO 4: PRESSURE TO DEPTH [cite: 2141-2150]
    ! =========================================================================
    FUNCTION ALGO4_DEPTH(P, LAT) RESULT(Z)
        REAL, INTENT(IN) :: P, LAT
        REAL :: Z, X, GR
        X = SIN(LAT / 57.29578)**2
        GR = 9.780318 * (1.0 + (5.2788E-3 + 2.36E-5*X)*X) + 1.092E-6*P
        Z = (((-1.82E-15*P + 2.279E-10)*P - 2.2512E-5)*P + 9.72659)*P
        Z = Z / GR
    END FUNCTION ALGO4_DEPTH

    ! =========================================================================
    ! [cite_start]ALGO 5: FREEZING POINT [cite: 2224-2236]
    ! =========================================================================
    FUNCTION ALGO5_FREEZE(S, P) RESULT(Tf)
        REAL, INTENT(IN) :: S, P
        REAL :: Tf
        Tf = (-0.0575*S) + (1.710523E-3*S*SQRT(S)) + (-2.154996E-4*S**2) + (-7.53E-4*P)
    END FUNCTION ALGO5_FREEZE

    ! =========================================================================
    ! [cite_start]ALGO 6: SPECIFIC HEAT (Cp) [cite: 2273-2332]
    ! =========================================================================
    FUNCTION ALGO6_SPEC_HEAT(S, T, P) RESULT(Cp)
        REAL, INTENT(IN) :: S, T, P
        REAL :: Cp, P_bar, SR
        REAL :: Cp_0, Del_Cp1, Del_Cp2, A, B, C
        
        P_bar = P / 10.0 ! Convert to Bars
        SR = SQRT(S)

        ! Cp(S,t,0) at Atmospheric Pressure
        Cp_0 = 4217.4 - 3.720283*T + 0.1412855*T**2 - 2.654387E-3*T**3 + 2.093236E-5*T**4 + &
               (-7.643575 + 0.1072763*T - 1.38385E-3*T**2)*S + &
               (0.1770383 - 4.07718E-3*T + 5.148E-5*T**2)*S*SR

        ! Pressure Term for Pure Water (Del_Cp1)
        A = -0.49592 + 1.45747E-2*T - 3.13885E-4*T**2 + 2.0357E-6*T**3 + 1.7168E-8*T**4
        B = 2.4931E-4 - 1.08645E-5*T + 2.87533E-7*T**2 - 4.0027E-9*T**3 + 2.2956E-11*T**4
        C = -5.422E-8 + 2.6380E-9*T - 6.5637E-11*T**2 + 6.136E-13*T**3
        Del_Cp1 = A*P_bar + B*P_bar**2 + C*P_bar**3

        ! Pressure Term for Salinity (Del_Cp2)
        A = (4.9247E-3 - 1.28315E-4*T + 9.802E-7*T**2 + 2.5941E-8*T**3 - 2.9179E-10*T**4) * S + &
            (-1.2331E-4 - 1.517E-6*T + 3.122E-8*T**2) * S*SR
        B = (-2.9558E-6 + 1.17054E-7*T - 2.3905E-9*T**2 + 1.8448E-11*T**3) * S + &
            (9.971E-8) * S*SR
        C = (5.540E-10 - 1.7682E-11*T + 3.513E-13*T**2) * S + &
            (-1.4300E-12 * T) * S*SR
        Del_Cp2 = A*P_bar + B*P_bar**2 + C*P_bar**3

        Cp = Cp_0 + Del_Cp1 + Del_Cp2
    END FUNCTION ALGO6_SPEC_HEAT

    ! =========================================================================
    ! [cite_start]ALGO 7: ADIABATIC LAPSE RATE (Gamma) [cite: 2455-2475]
    ! =========================================================================
    FUNCTION ALGO7_ADIABATIC(S, T, P) RESULT(Gamma)
        REAL, INTENT(IN) :: S, T, P
        REAL :: Gamma, DS
        DS = S - 35.0
        
        Gamma = 3.5803E-5 + 8.5258E-6*T - 6.8360E-8*T**2 + 6.6228E-10*T**3 + &
                (1.8932E-6 - 4.2393E-8*T)*DS + &
                (1.8741E-8 - 6.7795E-10*T + 8.7330E-12*T**2 - 5.4481E-14*T**3)*P + &
                (-1.1351E-10 + 2.7759E-12*T)*DS*P + &
                (-4.6206E-13 + 1.8676E-14*T - 2.1687E-16*T**2)*P**2
    END FUNCTION ALGO7_ADIABATIC

    ! =========================================================================
    ! [cite_start]ALGO 8: POTENTIAL TEMPERATURE (Theta) [cite: 2537-2562]
    ! Uses Runge-Kutta 4th Order Integration (Bryden 1973 / Fofonoff 1977)
    ! =========================================================================
    FUNCTION ALGO8_POT_TEMP(S, T_in, P_in, P_ref) RESULT(Theta)
        REAL, INTENT(IN) :: S, T_in, P_in, P_ref
        REAL :: Theta, H, Q1, Q2, Q3, Q4, T_local
        
        ! Step size (Delta P) = P_ref - P_in
        H = P_ref - P_in
        T_local = T_in

        ! RK4 Step 1
        Q1 = ALGO7_ADIABATIC(S, T_local, P_in)
        
        ! RK4 Step 2
        Q2 = ALGO7_ADIABATIC(S, T_local + 0.5*H*Q1, P_in + 0.5*H)
        
        ! RK4 Step 3
        Q3 = ALGO7_ADIABATIC(S, T_local + 0.5*H*Q2, P_in + 0.5*H)
        
        ! RK4 Step 4
        Q4 = ALGO7_ADIABATIC(S, T_local + H*Q3, P_in + H)

        ! Final Integration
        Theta = T_local + (H/6.0) * (Q1 + 2.0*Q2 + 2.0*Q3 + Q4)
    END FUNCTION ALGO8_POT_TEMP

    ! =========================================================================
    ! [cite_start]ALGO 9: SOUND SPEED (Chen & Millero) [cite: 2617-2665]
    ! =========================================================================
    FUNCTION ALGO9_SOUND(S, T, P) RESULT(U)
        REAL, INTENT(IN) :: S, T, P
        REAL :: U, Pb, SR, Cw, A, B, D
        Pb = P / 10.0
        SR = SQRT(ABS(S))
        Cw = 1402.388 + 5.03711*T - 0.0580852*T**2 + 3.3420E-4*T**3 - &
             1.47800E-6*T**4 + 3.1464E-9*T**5 + &
             (0.153563 + 6.8982E-4*T - 8.1788E-6*T**2 + 1.3621E-7*T**3 - 6.1185E-10*T**4)*Pb + &
             (3.1260E-5 - 1.7107E-6*T + 2.5974E-8*T**2 - 2.5335E-10*T**3 + 1.0405E-12*T**4)*Pb**2 + &
             (-9.7729E-9 + 3.8504E-10*T - 2.3643E-12*T**2)*Pb**3

        A = 1.389 - 1.262E-2*T + 7.164E-5*T**2 + 2.006E-6*T**3 - 3.21E-8*T**4 + &
            (9.4742E-5 - 1.2580E-5*T - 6.4885E-8*T**2 + 1.0507E-8*T**3 - 2.0122E-10*T**4)*Pb + &
            (-3.9064E-7 + 9.1041E-9*T - 1.6002E-10*T**2 + 7.988E-12*T**3)*Pb**2 + &
            (1.100E-10 + 6.649E-12*T - 3.389E-13*T**2)*Pb**3

        B = -1.922E-2 - 4.42E-5*T + (7.3637E-5 + 1.7945E-7*T)*Pb
        D = 1.727E-3 - 7.9836E-6*Pb
        U = Cw + A*S + B*S*SR + D*S**2
    END FUNCTION ALGO9_SOUND

END MODULE UNESCO_MATH

! ==============================================================================
! MAIN PROGRAM (DRIVER)
! ==============================================================================
PROGRAM MAIN_UNESCO
    USE UNESCO_MATH
    IMPLICIT NONE

    ! Variabel Input & Output
    INTEGER :: CHOICE
    REAL :: S_in, T_in, P_in, R_in, LAT_in, P_ref
    REAL :: Result_Val, Result_2
    LOGICAL :: RUNNING

    RUNNING = .TRUE.

    ! --- LOOP UTAMA (Materi Minggu 6) ---
    DO WHILE (RUNNING)
        PRINT *, " "
        PRINT *, "=========================================================="
        PRINT *, "      TUGAS OS2103: UNESCO ALGORITHMS (PAPER 44)          "
        PRINT *, "      Nama: [ISI NAMA KAMU] | NIM: [ISI NIM KAMU]         "
        PRINT *, "=========================================================="
        PRINT *, "1. Conductivity -> Salinity (Algo 1)"
        PRINT *, "2. Salinity -> Conductivity (Algo 2 - Iterative)"
        PRINT *, "3. Density & Spec. Volume (Algo 3)"
        PRINT *, "4. Pressure -> Depth (Algo 4)"
        PRINT *, "5. Freezing Point (Algo 5)"
        PRINT *, "6. Specific Heat (Algo 6)"
        PRINT *, "7. Adiabatic Lapse Rate (Algo 7)"
        PRINT *, "8. Potential Temperature (Algo 8 - Runge-Kutta)"
        PRINT *, "9. Sound Speed (Algo 9)"
        PRINT *, "0. EXIT"
        PRINT *, "=========================================================="
        WRITE (*, '(A)', ADVANCE='NO') " Masukkan Pilihan (0-9): "
        READ (*, *) CHOICE

        ! --- LOGIKA IF/ELSE (Materi Minggu 4) ---
        IF (CHOICE .EQ. 0) THEN
            RUNNING = .FALSE.
            PRINT *, "Program Selesai."
        
        ELSE IF (CHOICE .EQ. 1) THEN
            ! Input
            PRINT *, "Masukkan: Ratio (R), Temp (C), Pressure (db)"
            READ (*, *) R_in, T_in, P_in
            ! Hitung
            Result_Val = ALGO1_SAL78(R_in, T_in, P_in)
            ! Output (Sesuai Screenshot User)
            PRINT *, "Data ke-           1 :"
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   R (ratio)    T (C)       P (dbar)    Salinitas (PSU)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 100) R_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 2) THEN
            PRINT *, "Masukkan: Salinity (PSU), Temp (C), Pressure (db)"
            READ (*, *) S_in, T_in, P_in
            Result_Val = ALGO2_COND(S_in, T_in, P_in)
            PRINT *, "Data ke-           1 :"
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    T (C)       P (dbar)    R (ratio)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 100) S_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 3) THEN
            PRINT *, "Masukkan: Salinity (PSU), Temp (C), Pressure (db)"
            READ (*, *) S_in, T_in, P_in
            Result_Val = ALGO3_DENSITY(S_in, T_in, P_in) ! Density
            ! Spec Vol Anomaly approx (Optional, hanya density yg di output sini)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    T (C)       P (dbar)    Density (kg/m3)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 100) S_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 4) THEN
            PRINT *, "Masukkan: Pressure (db), Latitude (deg)"
            READ (*, *) P_in, LAT_in
            Result_Val = ALGO4_DEPTH(P_in, LAT_in)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Press (db)   Lat (deg)               Depth (m)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 101) P_in, LAT_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 5) THEN
            PRINT *, "Masukkan: Salinity (PSU), Pressure (db)"
            READ (*, *) S_in, P_in
            Result_Val = ALGO5_FREEZE(S_in, P_in)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    Press (db)              Freezing Pt (C)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 101) S_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 6) THEN
            PRINT *, "Masukkan: Salinity (PSU), Temp (C), Pressure (db)"
            READ (*, *) S_in, T_in, P_in
            Result_Val = ALGO6_SPEC_HEAT(S_in, T_in, P_in)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    T (C)       P (dbar)    Specific Heat (J/kg C)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 102) S_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 7) THEN
            PRINT *, "Masukkan: Salinity (PSU), Temp (C), Pressure (db)"
            READ (*, *) S_in, T_in, P_in
            Result_Val = ALGO7_ADIABATIC(S_in, T_in, P_in)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    T (C)       P (dbar)    Adiabatic Lapse (C/db)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 103) S_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 8) THEN
            PRINT *, "Masukkan: Salinity, Temp (in-situ), Pressure (in-situ), P_Ref"
            READ (*, *) S_in, T_in, P_in, P_ref
            Result_Val = ALGO8_POT_TEMP(S_in, T_in, P_in, P_ref)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   T (in-situ)  P (in-situ) P (ref)     Potential Temp (C)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 100) T_in, P_in, P_ref, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE IF (CHOICE .EQ. 9) THEN
            PRINT *, "Masukkan: Salinity, Temp, Pressure"
            READ (*, *) S_in, T_in, P_in
            Result_Val = ALGO9_SOUND(S_in, T_in, P_in)
            PRINT *, "----------------------------------------------------------------"
            PRINT *, "   Salinitas    T (C)       P (dbar)    Sound Spd (m/s)"
            PRINT *, "----------------------------------------------------------------"
            WRITE (*, 100) S_in, T_in, P_in, Result_Val
            PRINT *, "----------------------------------------------------------------"

        ELSE
            PRINT *, "Pilihan tidak valid."
        END IF

        ! --- FORMAT OUTPUT (Materi Minggu 7) ---
100     FORMAT (3X, F10.5, 2X, F10.5, 2X, F10.5, 4X, F10.5)
101     FORMAT (3X, F10.5, 2X, F10.5, 16X, F10.5)
102     FORMAT (3X, F10.5, 2X, F10.5, 2X, F10.5, 4X, F12.3) ! Utk Specific Heat (nilai ribuan)
103     FORMAT (3X, F10.5, 2X, F10.5, 2X, F10.5, 4X, E12.5) ! Utk Lapse Rate (nilai kecil)

        PRINT *, " "
    END DO

END PROGRAM MAIN_UNESCO