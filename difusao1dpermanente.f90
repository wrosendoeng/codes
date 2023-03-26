program difusao_1d_permanente

    implicit none
    
    ! Classificacao das variaveis
    character(len=1024)     :: metodo, nome_saida, str
    integer(4)              :: i, n, rc
    real(8)                 :: h, w, l, area, hconv, k, alfa, nu, t_int, t_ext, q0, sc, sp, psi
    real(8), dimension(:), allocatable :: x_p, dx_p, x_f, dx_f, ap, ae, aw, b, vetor_temperatura, &
    theta_numerico, theta_analitico, x_adim, e_abs, e_rel

    ! Dimensoes em mm -> Parelelepípedo de 1m x 1m x 0,1m
    h = 1.0d0       ! 1 metro
    w = 1.0d0       ! 1 metro
    l = 1.0d-1      ! 0,1 metro
    area = h * w    ! area transversal a conducao

    ! Definicao da condutividade termica de uma parede hipotética a T = 300 K (em W/(m.K))
    ! Fonte: Apendice A - livro "Transferencia de Calor - Incropera, 8a edicao"
    ! k = 1.0d0 ! tijolo refratario (tabela A.3)
    k = 6.05d1 ! aco carbono padrao (tabela A.1)
    ! k = 4.01d2 ! cobre padrao (tabela A.1) 

    ! Definicao do numero de Nusselt
    nu = 1.0d1 ! Definido pelo exercicio
    ! Definicao do Psi 
    psi = 5.0d0 ! Definido pelo exercicio
    hconv = nu*k/l ! h para conveccao
    ! Definicao das temperaturas do forno e do ambiente externo (em Kelvin)
    t_int = 1.273d3 ! 1000 graus Celsius
    t_ext = 3.0d2 ! 27 graus Celsius
    ! Definicao do fluxo inicial
    q0 = ((psi/l)**(2.0d0))*k*(t_ext-t_int)

    ! Definicao do fator de distribuicao da malha
    ! alfa > 1 --> concentracao dos pontos no início da malha
    ! alfa < 1 --> concentracao dos pontos no fim da malha
    alfa = 1.2d0
    ! numero de pontos a serem calculados
    n = 100 

    ! Equacao de conservacao
    ! d/dx(k*dT/dx) - q_ponto = 0
    ! q_ponto = q_ponto_inicial*(T - t_ext)/(t_int-t_ext)
    ! x = 0, T = t_int (condicao de contorno do tipo 01)
    ! x = L, -k*(dT/dx)|x=L = h*(T(L) - t_ext) (condicao de contorno do tipo 03)

    ! Definicao das posicoes dos pontos e das faces de acordo com o metodo A ou B
    metodo = "A"  ! tire o "!" da frente para usar o metodo A
    ! metodo = "B"  ! tire o "!" da frente para usar o metodo B

    select case (metodo)
    ! escolha pelo metodo A - meio volume em cada extremidade
        case("A")

        allocate(x_p(n), x_f(n+1), dx_p(n-1), dx_f(n), ap(n), ae(n), aw(n), b(n), &
        vetor_temperatura(n), theta_numerico(n), theta_analitico(n), x_adim(n), e_abs(n), e_rel(n))

        ! Calcular posicao dos pontos (x_p) e das faces (x_f):
        x_p(1) = 0.0d0
        do i = 2, n, 1
            x_p(i) = l*((i-1.0d0)/(n-1.0d0))**alfa
            x_f(i) = 0.5d0*(x_p(i) + x_p(i-1))
        end do
        x_f(1) = 0.0d0      ! 1a posicao da face
        x_f(n+1) = x_p(n)   ! N-esima posicao da face         
        ! Distancia entre os pontos
        do i = 1, n-1, 1
            dx_p(i) = x_p(i+1) - x_p(i)
            ! print *, dx_p(i)
        end do
        ! Distancia entre as faces
        do i = 1, n, 1
            dx_f(i) = x_f(i+1) - x_f(i)
            ! print *, dx_f(i)
        end do

    ! escolha pelo metodo B - volume zero nas extremidades
        case("B")
            
        allocate(x_p(n+1), x_f(n), dx_p(n), dx_f(n-1), ap(n), ae(n), aw(n), b(n), &
        vetor_temperatura(n), theta_numerico(n),theta_analitico(n), x_adim(n), e_abs(n), e_rel(n))
            
        ! Calcular posicao dos pontos (x_p) e das faces (x_f):
        x_p(1) = 0.0d0
        do i = 2, n, 1
            x_f(i-1) = l*((i-2.0d0)/(n-2.0d0))**alfa
            x_p(i) = 0.5d0*(x_f(i+1) + x_f(i))
        end do        
        ! Distancia entre os pontos
        do i = 1, n, 1
            dx_p(i) = x_p(i+1) - x_p(i)
        end do
        ! Distancia entre as faces
        do i = 1, n-1, 1
            dx_f(i) = x_f(i+1) - x_f(i)
        end do
    end select

    ! termos da fonte (ou sumidouro):
    ! Fonte => S = Sc + Sp*phi (em geral, Sp <= 0)
    ! Parte livre = Sc
    sc = q0*t_ext/(t_int-t_ext)
    ! Parte dependente = Sp (Garantindo que Sp seja negativo)
    sp = -q0/(t_int-t_ext)

    ! Vetores AP, AE, AW & B:
    ! vetor AP:
    ap(1) = 1.0d0                            ! necessario pela CC do tipo 1 (temperatura prescrita) 
    ap(n) = k/dx_p(n-1) + hconv - sp*dx_f(n) ! necessario pela CC do tipo 3 (conveccao e conducao) 
    ! vetor AE:
    ae(1) = 0.0d0 ! necessario pela CC do tipo 1 (temperatura prescrita) 
    ae(n) = 0.0d0 ! necessario pela CC do tipo 3 (conveccao e conducao) 
    ! vetor AW:
    aw(1) = 0.0d0 ! necessario pela CC do tipo 1 (temperatura prescrita) 
    aw(n) = k/dx_p(n-1)
    ! vetor B:
    b(1) = t_int                    ! necessario pela CC do tipo 1 (temperatura prescrita) 
    b(n) = sc*dx_f(n) + hconv*t_ext ! necessario pela CC do tipo 3 (conveccao e conducao) 
    
    ! iteracao para discretizar os termos dos coeficientes entre 
    do i = 2, n-1, 1
        ! VETOR AE
        ae(i) = k/dx_p(i)
        ! VETOR AW
        aw(i) = k/dx_p(i-1)
        ! VETOR AP
        ap(i) = ae(i) + aw(i) + sp*dx_f(i)
        ! VETOR B
        b(i) = -sc*dx_f(i)
    end do

    ! Variavel de interesse -> Temperatura
    vetor_temperatura = 0.0d0 ! inicializando vetor de resposta
    
    !! Chamando a subrotina do algoritmo TDMA (ou Algoritmo de Thomas)
    write(str,"(I3)") n
    write(nome_saida, '("trabalho01_psi5_alfa12_", (a), "pontos.txt")' ) trim(str)
    open(unit=44,file=nome_saida,action='write',status='replace')
    call solve_tridiag(ap,ae,aw,b,vetor_temperatura,n)

    ! adimensionalizacao das propriedades:
    ! X_adim -> Comprimento adimensional (X_adim = x/L)
    x_adim = x_p/l
    ! Theta -> temperatura adimensional (theta = (T - T_ext)/(T_int - T_ext)) obtido numericamente
    theta_numerico = (vetor_temperatura-t_ext)/(t_int-t_ext)
    ! Theta -> temperatura adimensional (theta = 1 + (Psi*tanh(Psi)-Nu)/(Nu*tanh(Psi)+Psi)*tanh(Psi))
    ! Calculando o termo analitico de acordo com a fonte adimensionalizada (Psi)
    if (psi == 0.0d0) then
        theta_analitico = 1.0d0 - nu*x_adim/(nu + 1.0d0) ! sem geracao ou extracao de calor
    else
        theta_analitico = 1.0d0 - (psi*tanh(psi)+nu)*tanh(psi*x_adim)/(nu*tanh(psi)+psi) ! com geracao
    end if
    
    ! calculo do erro absoluto e relativo:
    e_abs = abs(theta_numerico-theta_analitico)
    e_rel = e_abs/theta_analitico

    do i = 1, n, 1
        ! Na sequencia: posicao do ponto em x; temperatura; posicao adimensionalizada; theta numerico; theta analitico; erro absoluto e erro relativo
        write(44,'(*(f15.8))') x_p(i), vetor_temperatura(i), x_adim(i), theta_numerico(i), &
        theta_analitico(i), e_abs(i), e_rel(i)
    end do

    ! Dealocacao de variaveis
    deallocate(x_p, dx_p, x_f, dx_f, ap, ae, aw, b, vetor_temperatura, theta_numerico, e_abs, e_rel)

    ! Fechar a unidade de salvamento do arquivo
    close(44)

    ! Inserir as subrotinas e modulos
    contains

    subroutine solve_tridiag(a,b,c,d,phi,nponto)
        implicit none
    
        integer,intent(in)                      :: nponto
        real(8),dimension(nponto),intent(in)    :: a, b, c, d 
        real(8),dimension(nponto),intent(out)   :: phi
        real(8),dimension(nponto)               :: p, q
        integer                                 :: iter

    ! inicializar os termos P e Q do algoritmo TDMA:
        p(1) = b(1)/a(1)
        q(1) = d(1)/a(1)

    ! Resolver os termos P e Q
        do iter = 2, nponto, 1
            p(iter) = b(iter)/(a(iter) - c(iter)*p(iter-1))
            q(iter) = (d(iter) + c(iter)*q(iter-1))/(a(iter) - c(iter)*p(iter-1))
        end do

        phi(nponto) = q(nponto)

    ! Resolver o vetor resposta:
        do iter = nponto-1, 1, -1
            phi(iter) = q(iter) + p(iter)*phi(iter+1)
        end do

    end subroutine solve_tridiag
  
end program difusao_1d_permanente