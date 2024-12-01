module sbs

use scfs

contains

subroutine generate_strs(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, perm_of, perm_ob, & 
    ph_of, fac_of, fac_ob, id_f, phase, binom, num_th)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, ncf
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: conff(ncf), confb(ncf)
    integer(8), intent(in) :: perm_of(nof), perm_ob(nob), ph_of(nof)
    complex(8), intent(in) :: fac_of(nof), fac_ob(nob)
    integer(8), intent(out) :: id_f(ncf)
    complex(8), intent(out) :: phase(ncf)
    integer(8), intent(in) :: num_th 

    integer(8) :: i, of, of1, cff, cff1, cffp, sgn, cff0 
    integer(8) :: ob, ob1, cfb, nb(nob), nb1(nob)
    complex(8) :: phi

    cff0 = 0 
    do of = 0, nof - 1
        if (ph_of(of + 1) == 1) cff0 = ibset(cff0, of)
    end do
    call omp_set_num_threads(num_th)
    !$omp parallel shared(nof, nob, ncf, lid, rid, conff, confb, perm_of, perm_ob, ph_of, fac_of, fac_ob, id_f, phase, cff0), &
    !$omp& private(i, of, of1, cff, cff1, cffp, sgn, ob, ob1, cfb, nb, nb1, phi)
    !$omp do 
    do i = 1, ncf
        phi = 1.
        cfb = confb(i)
        nb = decode_nb(nob, nebm, cfb, binom)
        nb1 = 0 
        do ob = 1, nob 
            ob1 = perm_ob(ob)
            phi = phi * fac_ob(ob) ** nb(ob)
            nb1(ob1) = nb(ob)
        end do

        cff = conff(i)
        cff1 = cff0
        cffp = 0
        do of = nof - 1, 0, -1
            of1 = perm_of(of + 1) - 1 
            if (ibits(cff, of, 1_8) == 0) cycle
            phi = phi * fac_of(of + 1)
            cff1 = ieor(cff1, ibset(0_8, of1))
            cffp = ieor(cffp, ibits(cff1, 0_8, of1))
        end do
        sgn = 0
        do of = 0, nof - 1
            sgn = ieor(sgn, ibits(cffp, of, 1_8))
        end do
        if (sgn == 1) phi = -phi
        id_f(i) = search_scf(nof, nob, norf, norb, nebm, lid, rid, cff1, nb1, binom)
        phase(i) = phi
    end do
    !$omp end do
    !$omp end parallel
end subroutine

subroutine generate_sbs_cfgr(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, nqnz, qnz_s, cyc, &
        perm_of, perm_ob, ph_of, fac_of, fac_ob, szz, dim, cfgr, cffac, binom, num_th, disp_std)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb, nebm, ncf
    integer(8), intent(in) :: binom(nebm + 1, nob + 1)
    integer(8), intent(in) :: lid(ishft(binom(nebm + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid(ishft(binom(nebm + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: conff(ncf), confb(ncf)
    integer(8), intent(in) :: nqnz, cyc(nqnz), szz
    complex(8), intent(in) :: qnz_s(nqnz)
    integer(8), intent(in) :: perm_of(nof, nqnz), perm_ob(nob, nqnz), ph_of(nof, nqnz)
    complex(8), intent(in) :: fac_of(nof, nqnz), fac_ob(nob, nqnz)
    integer(8), intent(out) :: dim, cfgr(ncf)
    complex(8), intent(out) :: cffac(ncf)
    integer(8) :: id(szz)
    complex(8) :: fac(szz)
    integer(8), allocatable :: transf(:,:)
    complex(8), allocatable :: phase(:,:)
    integer(8) :: i, j, k, sz
    complex(8) :: wt0
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    allocate(transf(ncf, nqnz))
    allocate(phase(ncf, nqnz))

    if (disp_std == 1) print *, 'Generating transformations'
    do j = 1, nqnz
        if (abs(qnz_s(j)) < 1.d-6) cycle
        call generate_strs(nof, nob, norf, norb, nebm, ncf, lid, rid, conff, confb, perm_of(:, j), perm_ob(:, j), ph_of(:, j), &
            fac_of(:, j), fac_ob(:, j), transf(:, j), phase(:, j), binom, num_th)
    end do
    if (disp_std == 1) print *, 'Generating basis cfgr'
    dim = 0
    cfgr = -1
    do i = 1, ncf
        id(1) = i 
        if (cfgr(id(1)) /= -1) cycle 
        fac(1) = 1.
        sz = 1
        do j = 1, nqnz 
            if (qnz_s(j) == 0) cycle 
            do k = 1, sz * (cyc(j) - 1)
                id(k + sz) = transf(id(k), j)
                fac(k + sz) = fac(k) * phase(id(k), j) * qnz_s(j)
            end do
            sz = sz * cyc(j)
        end do
        wt0 = 0
        do j = 1, sz 
            if (id(1) == id(j)) wt0 = wt0 + fac(j)
        end do 
        if (abs(wt0) < 1.d-6) cycle 
        dim = dim + 1
        cfgr(id(1 : sz)) = dim 
        cffac(id(1 : sz)) = fac(1 : sz) * sqrt(wt0 / sz)
    end do

    deallocate(transf)
    deallocate(phase)

end subroutine

subroutine generate_sbs_grel(ncf, szz, dim, cfgr, grel, grsz, disp_std)
    implicit none
    integer(8), intent(in) :: ncf, dim, cfgr(ncf), szz
    integer(8), intent(out) :: grel(szz, dim), grsz(dim)
    integer(8) :: i, g
    integer(8), intent(in) :: disp_std

    if (disp_std == 1) print *, 'Generating basis grel, dimension :', dim
    grel = -1 
    grsz = 0 
    do i = 1, ncf 
        g = cfgr(i)
        if (g == -1) cycle 
        grsz(g) = grsz(g) + 1
        grel(grsz(g), g) = i
    end do
    if (disp_std == 1) print *, 'Generating basis finish'
end subroutine

end module