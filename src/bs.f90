module bs

use cfs

contains

subroutine generate_trs(no, nor, ncf, lid, rid, conf, perm_o, ph_o, fac_o, id_f, phase, num_th)
    implicit none
    integer(8), intent(in) :: no, nor, ncf
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1), rid(ibset(0_8, nor)), conf(ncf)
    integer(8), intent(in) :: perm_o(no), ph_o(no)
    complex(8), intent(in) :: fac_o(no)
    integer(8), intent(out) :: id_f(ncf)
    complex(8), intent(out) :: phase(ncf)
    integer(8), intent(in) :: num_th 

    integer(8) :: i, o, o1, cf, cf1, cfp, ne, sgn, cf0 
    complex(8) :: phi

    cf0 = 0 
    do o = 0, no - 1
        if (ph_o(o + 1) == 1) cf0 = ibset(cf0, o)
    end do
    call omp_set_num_threads(num_th)
    !$omp parallel shared(no, nor, perm_o, fac_o, ncf, lid, rid, conf, id_f, phase, cf0), &
    !$omp& private(i, o, o1, cf, cf1, cfp, ne, sgn, phi)
    !$omp do 
    do i = 1, ncf
        cf = conf(i)
        cf1 = cf0
        ne = 0
        phi = 1.
        cfp = 0
        do o = no - 1, 0, -1
            o1 = perm_o(o + 1) - 1 
            if (ibits(cf, o, 1_8) == 0) cycle
            phi = phi * fac_o(o + 1)
            cf1 = ieor(cf1, ibset(0_8, o1))
            cfp = ieor(cfp, ibits(cf1, 0_8, o1))
        end do
        sgn = 0
        do o = 0, no - 1
            sgn = ieor(sgn, ibits(cfp, o, 1_8))
        end do
        if (sgn == 1) phi = -phi
        id_f(i) = search_conf(no, nor, lid, rid, cf1)
        phase(i) = phi
    end do
    !$omp end do
    !$omp end parallel
end subroutine

subroutine generate_bs_cfgr(no, nor, ncf, lid, rid, conf, nqnz, qnz_s, cyc, perm_o, ph_o, fac_o, szz, dim, cfgr, cffac, &
        num_th, silence_std)
    implicit none
    integer(8), intent(in) :: no, nor, ncf
    integer(8), intent(in) :: lid(ibset(0_8, no - nor) + 1), rid(ibset(0_8, nor)), conf(ncf)
    integer(8), intent(in) :: nqnz, cyc(nqnz), szz
    complex(8), intent(in) :: qnz_s(nqnz)
    integer(8), intent(in) :: perm_o(no, nqnz), ph_o(no, nqnz)
    complex(8), intent(in) :: fac_o(no, nqnz)
    integer(8), intent(out) :: dim, cfgr(ncf)
    complex(8), intent(out) :: cffac(ncf)
    integer(8) :: id(szz)
    complex(8) :: fac(szz)
    integer(8), allocatable :: transf(:,:)
    complex(8), allocatable :: phase(:,:)
    integer(8) :: i, j, k, sz
    real(8) :: wt0
    integer(8), intent(in) :: num_th 
    logical, intent(in) :: silence_std

    allocate(transf(ncf, nqnz))
    allocate(phase(ncf, nqnz))

    if (.not. silence_std) print *, 'Generating transformations'
    do j = 1, nqnz
        if (abs(qnz_s(j)) < 1.d-6) cycle
        call generate_trs(no, nor, ncf, lid, rid, conf, perm_o(:, j), ph_o(:, j), fac_o(:, j), transf(:, j), phase(:, j), num_th)
    end do
    if (.not. silence_std) print *, 'Generating basis cfgr'
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

subroutine generate_bs_grel(ncf, szz, dim, cfgr, grel, grsz, silence_std)
    implicit none
    integer(8), intent(in) :: ncf, dim, cfgr(ncf), szz
    integer(8), intent(out) :: grel(szz, dim), grsz(dim)
    integer(8) :: i, g
    logical, intent(in) :: silence_std

    if (.not. silence_std) print *, 'Generating basis grel'
    grel = -1 
    grsz = 0 
    do i = 1, ncf 
        g = cfgr(i)
        if (g == -1) cycle 
        grsz(g) = grsz(g) + 1
        grel(grsz(g), g) = i
    end do
    if (.not. silence_std) print *, 'Generating basis finish'
end subroutine

end module