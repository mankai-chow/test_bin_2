module scfs 
    
contains 

subroutine count_scfs(nof, nob, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: nof, nob, nqnu
    integer(8), intent(in) :: qnu_s(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(out) :: lid(ibset(0_8, nof) + 1)
    integer(8), intent(out) :: ncf
    integer(8) :: qnu_1(nqnu), i, j
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Counting configurations start, number of threads :', omp_get_max_threads()
    !$omp parallel shared(nof, nob, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, disp_std) private(qnu_1, i, j)
    !$omp do
    do i = 0, ibset(0_8, nof) - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Counting configurations', &
                i + 1, '*', omp_get_max_threads(), '/', ibset(0_8, nof)
        end if 
        qnu_1 = qnu_s
        do j = 1, nof
            if (ibits(i, j - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_of(j, :)
        end do
        lid(i + 2) = 0
        call count_scfs_rec(nob, 1_8, nqnu, qnu_1, qnu_ob, modul, lid(i + 2))
    end do
    !$omp end do
    !$omp end parallel
    lid(1) = 0
    do i = 1, ibset(0_8, nof)
        lid(i + 1) = lid(i + 1) + lid(i)
    end do
    ncf = lid(ibset(0_8, nof) + 1)
    if (disp_std == 1) print *, 'Counting configurations finish, total number :', ncf
    
end subroutine

recursive subroutine count_scfs_rec(nob, noc, nqnu, qnu_1, qnu_ob, modul, ct)
    implicit none
    integer(8), intent(in) :: nob, noc, nqnu
    integer(8), intent(in) :: qnu_1(nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(inout) :: ct
    integer(8) :: i, cyc1, qnu_c(nqnu)

    if (noc == nob + 1) then 
        do i = 1, nqnu
            cyc1 = modul(i)
            if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) /= 0) return
            if (cyc1 == 1 .and. qnu_1(i) /= 0) return
        end do
        ct = ct + 1 
        return
    end if
    
    qnu_c = qnu_1
    do 
        call count_scfs_rec(nob, noc + 1, nqnu, qnu_c, qnu_ob, modul, ct)
        qnu_c = qnu_c - qnu_ob(noc, :)
        do i = 1, nqnu 
            if (modul(i) == 1 .and. qnu_c(i) < 0) return 
        end do
    end do 
end subroutine 

subroutine generate_scfs(nof, nob, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, rid, conff, confb, num_th, disp_std)
    use omp_lib
    implicit none
    integer(8), intent(in) :: nof, nob, nebm, nqnu, ncf
    integer(8), intent(in) :: qnu_s(nqnu), qnu_of(nof, nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(ibset(0_8, nof) + 1)
    integer(8), intent(out) :: conff(ncf), confb(ncf), rid(ibset(0_8, nob + nebm))
    integer(8) :: qnu_1(nqnu), i, j, ct, tmp(nob)
    integer(8), intent(in) :: num_th 
    integer(8), intent(in) :: disp_std

    call omp_set_num_threads(num_th)
    if (disp_std == 1) print *, 'Generating configurations start'
    rid = 0
    !$omp parallel shared(nof, nob, nebm, nqnu, qnu_s, qnu_of, qnu_ob, modul, ncf, lid, rid, conff, confb, disp_std), &
    !$omp& private(qnu_1, i, j, ct, tmp)
    !$omp do
    do i = 0, ibset(0_8, nof) - 1
        if (mod(i + 1, 10000) == 0 .and. disp_std == 1) then 
            if (omp_get_thread_num() == 0) print *, 'Generating configurations', &
                i + 1, '*', omp_get_max_threads(), '/', ibset(0_8, nof)
        end if 
        if (lid(i + 1) == lid(i + 2)) cycle 
        qnu_1 = qnu_s
        do j = 1, nof
            if (ibits(i, j - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_of(j, :)
        end do
        ct = lid(i + 1)
        tmp = 0
        call generate_scfs_rec(nof, nob, nebm, 1_8, nqnu, qnu_1, qnu_ob, modul, ncf, ct, i, tmp, lid, rid, conff, confb)
    end do
    !$omp end do
    !$omp end parallel
    if (disp_std == 1)  print *, 'Generating configurations finish'
end subroutine

recursive subroutine generate_scfs_rec(nof, nob, nebm, noc, nqnu, qnu_1, qnu_ob, modul, ncf, ct, cff, tmp, lid, rid, conff, confb)
    implicit none
    integer(8), intent(in) :: nof, nob, nebm, noc, nqnu, tmp(nob), ncf, cff
    integer(8), intent(in) :: qnu_1(nqnu), qnu_ob(nob, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(ibset(0_8, nof) + 1)
    integer(8), intent(inout) :: ct, conff(ncf), confb(ncf), rid(ibset(0_8, nob + nebm))
    integer(8) :: i, cyc1, qnu_c(nqnu), tmp_1(nob), cfb

    if (noc == nob + 1) then 
        do i = 1, nqnu
            cyc1 = modul(i)
            if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) /= 0) return
            if (cyc1 == 1 .and. qnu_1(i) /= 0) return
        end do
        ct = ct + 1 
        conff(ct) = cff 
        cfb = encode_nb(nob, tmp)
        confb(ct) = cfb
        if (rid(cfb + 1) == 0) then
            rid(cfb + 1) = ct - lid(cff + 1)
        end if
        return
    end if

    qnu_c = qnu_1 
    tmp_1 = tmp
    do 
        call generate_scfs_rec(nof, nob, nebm, noc + 1, nqnu, qnu_c, qnu_ob, modul, ncf, ct, cff, tmp_1, lid, rid, conff, confb)
        tmp_1(noc) = tmp_1(noc) + 1
        qnu_c = qnu_c - qnu_ob(noc, :)
        do i = 1, nqnu 
            if (modul(i) == 1 .and. qnu_c(i) < 0) return 
        end do        
    end do
end subroutine

function encode_nb(nob, nb) result(cfb)
    implicit none
    integer(8), intent(in) :: nob, nb(nob)
    integer(8) :: cfb 
    integer(8) :: i, pos 
    pos = -1
    cfb = 0
    do i = 1, nob 
        pos = pos + 1 + nb(i)
        cfb = ibset(cfb, pos)
    end do
end function 

function decode_nb(nob, cfb) result(nb)
    implicit none
    integer(8), intent(in) :: nob, cfb 
    integer(8) :: nb(nob)
    integer(8) :: i, pos, ob 
    pos = -1
    ob = 1
    do 
        pos = pos + 1
        if (ibits(cfb, pos, 1) == 0) cycle 
        nb(ob) = pos 
        ob = ob + 1
        if (ob == nob + 1) exit
    end do
    do i = nob, 2, -1
        nb(i) = nb(i) - nb(i - 1) - 1
    end do
end function 

function search_scf(nof, nob, lid, rid, cff, cfb) result(i)
    implicit none
    integer(8), intent(in) :: nof, nob
    integer(8), intent(in) :: lid(ibset(0_8, nof) + 1), rid(ibset(0_8, nob))
    integer(8), intent(in) :: cff, cfb
    integer(8):: i
    i = lid(cff + 1) + rid(cfb + 1)
end function

end