program teste

    implicit none

    real(8) :: numero, valor
    real(8), dimension(3) :: vetor, vetor2, res

    vetor = (/3.0d0, 4.0d0, 5.0d0/)
    vetor2 = (/1.5d0, -2.6d0, -3.7d0/)
    
    res = abs(vetor-vetor2)/vetor2

    print *, res

end program teste