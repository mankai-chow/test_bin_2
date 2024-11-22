program main 
    use diag_re
    implicit none
    
    integer(8), parameter :: dim = 20
    integer(8), parameter :: sym_q = 0
    integer(8), parameter :: nel = 40
    integer(8), parameter :: nst = 5
    real(8), parameter :: tol = 1.0d-8
    integer(8), parameter :: ncv_in = 15
    integer(8), parameter :: num_th = 8
    integer(8), parameter :: display_std = 1

    integer(8) :: i
    integer(8) :: colptr(dim + 1), rowid(nel)
    real(8) :: elval(nel)
    real(8) :: eigval(nst + 1), eigvec(dim, nst + 1)

    do i = 1, dim 
        colptr(i) = 2 * i - 1
        rowid(i * 2 - 1) = i - 1
        rowid(i * 2) = i + 1
        if (i == 1) rowid(i * 2 - 1) = dim 
        if (i == dim) rowid(i * 2) = 1
    end do 
    colptr(dim + 1) = nel + 1
    elval = 1.0

    call diagonalisation_re(dim, sym_q, nel, colptr, rowid, elval, nst, tol, ncv_in, eigval, eigvec, num_th, display_std)

    print '(1f16.8)', eigval(:)

end program