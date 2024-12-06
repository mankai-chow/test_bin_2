module sent 

use scfs 

contains

subroutine decomp_sbasis(nof, nob, norf, norb, &
    nebm_0, ncf_0, dim_0, conff_0, confb_0, lid_0, rid_0, szz_0, cfgr_0, cffac_0, grel_0, grsz_0, binom_0, &
    nebm_a, ncf_a, dim_a, conff_a, confb_a, lid_a, rid_a, szz_a, cfgr_a, cffac_a, grel_a, grsz_a, binom_a, &
    nebm_b, ncf_b, dim_b, conff_b, confb_b, lid_b, rid_b, szz_b, cfgr_b, cffac_b, grel_b, grsz_b, binom_b, &
    amp_ofa, amp_oba, amp_ofb, amp_obb, binom_c, st_0, st_dcp)
    implicit none
    integer(8), intent(in) :: nof, nob, norf, norb

    integer(8), intent(in) :: nebm_0, ncf_0, dim_0, szz_0
    integer(8), intent(in) :: conff_0(ncf_0), confb_0(ncf_0)
    integer(8), intent(in) :: binom_0(nebm_0 + 1, nob + 1)
    integer(8), intent(in) :: lid_0(ishft(binom_0(nebm_0 + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid_0(ishft(binom_0(nebm_0 + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: cfgr_0(ncf_0), grel_0(szz_0, dim_0), grsz_0(dim_0)
    complex(8), intent(in) :: cffac_0(ncf_0)

    integer(8), intent(in) :: nebm_a, ncf_a, dim_a, szz_a
    integer(8), intent(in) :: conff_a(ncf_a), confb_a(ncf_a)
    integer(8), intent(in) :: binom_a(nebm_a + 1, nob + 1)
    integer(8), intent(in) :: lid_a(ishft(binom_a(nebm_a + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid_a(ishft(binom_a(nebm_a + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: cfgr_a(ncf_a), grel_a(szz_a, dim_a), grsz_a(dim_a)
    complex(8), intent(in) :: cffac_a(ncf_a)

    integer(8), intent(in) :: nebm_b, ncf_b, dim_b, szz_b
    integer(8), intent(in) :: conff_b(ncf_b), confb_b(ncf_b)
    integer(8), intent(in) :: binom_b(nebm_b + 1, nob + 1)
    integer(8), intent(in) :: lid_b(ishft(binom_b(nebm_b + 1, nob - norb + 1) + 1, nof - norf) + 1)
    integer(8), intent(in) :: rid_b(ishft(binom_b(nebm_b + 1, norb + 1) + 1, norf) + 1)
    integer(8), intent(in) :: cfgr_b(ncf_b), grel_b(szz_b, dim_b), grsz_b(dim_b)
    complex(8), intent(in) :: cffac_b(ncf_b)

    complex(8), intent(in) :: amp_ofa(nof), amp_oba(nob), amp_ofb(nof), amp_obb(nob), st_0(dim_0)
    integer(8), intent(in) :: binom_c(nebm_a + 1, nebm_b + 1)
    complex(8), intent(out) :: st_dcp(dim_b, dim_a)

    integer(8) :: ga, gb, g0, ea, eb, ia, ib, i0, o, ph
    integer(8) :: cffa, cffb, cff0, cffp, cfb0, nba(nob), nbb(nob), nb0(nob)
    complex(8) :: amp_cfa(ncf_a), amp_cfb(ncf_b), fac2

    do ia = 1, ncf_a 
        amp_cfa(ia) = 1 
        cffa = conff_a(ia)
        do o = 0, nof - 1
            if (ibits(cffa, o, 1_8) == 1_8) amp_cfa(ia) = amp_cfa(ia) * amp_ofa(o + 1)
        end do
        nba = decode_nb(nob, nebm_a, confb_a(ia), binom_a)
        do o = 1, nob 
            amp_cfa(ia) = amp_cfa(ia) * (amp_oba(o) ** nba(o))
        end do
    end do
    do ib = 1, ncf_b 
        amp_cfb(ib) = 1 
        cffb = conff_b(ib)
        do o = 0, nof - 1
            if (ibits(cffb, o, 1_8) == 1_8) amp_cfb(ib) = amp_cfb(ib) * amp_ofb(o + 1)
        end do
        nbb = decode_nb(nob, nebm_b, confb_b(ib), binom_b)
        do o = 1, nob 
            amp_cfa(ib) = amp_cfa(ib) * (amp_oba(o) ** nba(o))
        end do
    end do

    do ga = 1, dim_a 
        do gb = 1, dim_b 
            st_dcp(gb, ga) = 0
            do ea = 1, grsz_a(ga) 
                ia = grel_a(ea, ga)
                cffa = conff_a(ia)
                nba = decode_nb(nob, nebm_a, confb_a(ia), binom_a)
                do eb = 1, grsz_b(gb) 
                    ib = grel_b(eb, gb) 
                    cffb = conff_b(ib)
                    nbb = decode_nb(nob, nebm_b, confb_b(ib), binom_b)
                    if (iand(cffa, cffb) /= 0) cycle
                    cff0 = ior(cffa, cffb) 
                    nb0 = nba + nbb
                    i0 = search_scf(nof, nob, norf, norb, nebm_0, lid_0, rid_0, cff0, nb0, binom_0)
                    if (i0 <= 0 .or. i0 > ncf_0) cycle
                    if (cff0 /= conff_0(i0)) cycle
                    cfb0 = encode_nb(nob, nebm_0, nb0, binom_0)
                    if (cfb0 /= confb_0(i0)) cycle
                    g0 = cfgr_0(i0)
                    if (g0 == -1) cycle
                    cffp = 0
                    do o = 0, nof - 1
                        if (ibits(cffa, o, 1_8) == 0_8) cycle
                        cffp = ieor(cffp, ibits(cffb, 0_8, o))
                    end do
                    ph = 1
                    do o = 0, nof - 1
                        if (ibits(cffp, o, 1_8) == 0_8) ph = -ph
                    end do
                    fac2 = 1.0
                    do o = 1, nob 
                        fac2 = fac2 * binom_c(nba(o) + 1, nbb(o) + 1)
                    end do
                    st_dcp(gb, ga) = st_dcp(gb, ga) + st_0(g0) * ph * sqrt(fac2) * &
                        cffac_0(i0) * cffac_a(ia) * cffac_b(ib) * amp_cfa(ia) * amp_cfb(ib)
                end do
            end do
        end do
    end do
end subroutine

end module