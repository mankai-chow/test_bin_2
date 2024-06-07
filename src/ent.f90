module ent 

use cfs 

contains

subroutine decomp_basis(no, nor, &
    ncf_0, dim_0, conf_0, lid_0, rid_0, szz_0, cfgr_0, cffac_0, grel_0, grsz_0, &
    ncf_a, dim_a, conf_a, lid_a, rid_a, szz_a, cfgr_a, cffac_a, grel_a, grsz_a, &
    ncf_b, dim_b, conf_b, lid_b, rid_b, szz_b, cfgr_b, cffac_b, grel_b, grsz_b, &
    amp_oa, amp_ob, st_0, st_dcp)
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_0, dim_0, szz_0
    integer(8), intent(in) :: conf_0(ncf_0), lid_0(ibset(0_8, no - nor)), rid_0(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_0(ncf_0), grel_0(szz_0, dim_0), grsz_0(dim_0)
    complex(8), intent(in) :: cffac_0(ncf_0)

    integer(8), intent(in) :: ncf_a, dim_a, szz_a
    integer(8), intent(in) :: conf_a(ncf_a), lid_a(ibset(0_8, no - nor)), rid_a(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_a(ncf_a), grel_a(szz_a, dim_a), grsz_a(dim_a)
    complex(8), intent(in) :: cffac_a(ncf_a)

    integer(8), intent(in) :: ncf_b, dim_b, szz_b
    integer(8), intent(in) :: conf_b(ncf_b), lid_b(ibset(0_8, no - nor) + 1), rid_b(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_b(ncf_b), grel_b(szz_b, dim_b), grsz_b(dim_b)
    complex(8), intent(in) :: cffac_b(ncf_b)

    complex(8), intent(in) :: amp_oa(no), amp_ob(no), st_0(dim_0)
    complex(8), intent(out) :: st_dcp(dim_b, dim_a)

    integer(8) :: ga, gb, g0, ea, eb, ia, ib, i0, cfa, cfb, cf0, cfp, o, ph
    complex(8) :: fac, amp_cfa(ncf_a), amp_cfb(ncf_b)

    do ia = 1, ncf_a 
        amp_cfa(ia) = 1 
        cfa = conf_a(ia)
        do o = 0, no - 1
            if (ibits(cfa, o, 1_8) == 1_8) amp_cfa(ia) = amp_cfa(ia) * amp_oa(o + 1)
        end do
    end do
    do ib = 1, ncf_b 
        amp_cfb(ib) = 1 
        cfb = conf_b(ib)
        do o = 0, no - 1
            if (ibits(cfb, o, 1_8) == 1_8) amp_cfb(ib) = amp_cfb(ib) * amp_ob(o + 1)
        end do
    end do

    do ga = 1, dim_a 
        do gb = 1, dim_b 
            st_dcp(gb, ga) = 0
            do ea = 1, grsz_a(ga) 
                ia = grel_a(ea, ga)
                cfa = conf_a(ia)
                do eb = 1, grsz_b(gb) 
                    ib = grel_b(eb, gb) 
                    cfb = conf_b(ib)
                    if (iand(cfa, cfb) /= 0) cycle
                    cf0 = ior(cfa, cfb) 
                    i0 = search_conf(no, nor, lid_0, rid_0, cf0)
                    if (i0 <= 0 .or. i0 > ncf_0) cycle
                    if (cf0 /= conf_0(i0)) cycle
                    g0 = cfgr_0(i0)
                    if (g0 == -1) cycle
                    cfp = 0
                    do o = 0, no - 1
                        if (ibits(cfa, o, 1_8) == 0_8) cycle
                        cfp = ieor(cfp, ibits(cfb, 0_8, o))
                    end do
                    ph = 1
                    do o = 0, no - 1
                        if (ibits(cfp, o, 1_8) == 0_8) ph = -ph
                    end do
                    st_dcp(gb, ga) = st_dcp(gb, ga) + &
                        st_0(g0) * ph * cffac_0(i0) * cffac_a(ia) * cffac_b(ib) * amp_cfa(ia) * amp_cfb(ib)
                end do
            end do
        end do
    end do
end subroutine

subroutine combine_basis(no, nor, &
    ncf_0, dim_0, conf_0, lid_0, rid_0, szz_0, cfgr_0, cffac_0, grel_0, grsz_0, &
    ncf_a, dim_a, conf_a, lid_a, rid_a, szz_a, cfgr_a, cffac_a, grel_a, grsz_a, &
    ncf_b, dim_b, conf_b, lid_b, rid_b, szz_b, cfgr_b, cffac_b, grel_b, grsz_b, &
    st_dcp, st_0)
    implicit none
    integer(8), intent(in) :: no, nor 

    integer(8), intent(in) :: ncf_0, dim_0, szz_0
    integer(8), intent(in) :: conf_0(ncf_0), lid_0(ibset(0_8, no - nor)), rid_0(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_0(ncf_0), grel_0(szz_0, dim_0), grsz_0(dim_0)
    complex(8), intent(in) :: cffac_0(ncf_0)

    integer(8), intent(in) :: ncf_a, dim_a, szz_a
    integer(8), intent(in) :: conf_a(ncf_a), lid_a(ibset(0_8, no - nor)), rid_a(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_a(ncf_a), grel_a(szz_a, dim_a), grsz_a(dim_a)
    complex(8), intent(in) :: cffac_a(ncf_a)

    integer(8), intent(in) :: ncf_b, dim_b, szz_b
    integer(8), intent(in) :: conf_b(ncf_b), lid_b(ibset(0_8, no - nor) + 1), rid_b(ibset(0_8, nor))
    integer(8), intent(in) :: cfgr_b(ncf_b), grel_b(szz_b, dim_b), grsz_b(dim_b)
    complex(8), intent(in) :: cffac_b(ncf_b)

    complex(8), intent(in) :: st_dcp(dim_b, dim_a)
    complex(8), intent(out) :: st_0(dim_0)

    integer(8) :: ga, gb, g0, ea, eb, ia, ib, i0, cfa, cfb, cf0, cfp, o, ph
    complex(8) :: fac

    print *, 'WARNING : THIS FUNCTION MAY HAVE TROUBLE DEALING WITH DOUBLE OCCUPANCY'
    st_0 = 0
    do ga = 1, dim_a 
        do gb = 1, dim_b 
            do ea = 1, grsz_a(ga) 
                ia = grel_a(ea, ga)
                cfa = conf_a(ia)
                do eb = 1, grsz_b(gb) 
                    ib = grel_b(eb, gb) 
                    cfb = conf_b(ib)
                    if (iand(cfa, cfb) /= 0) cycle
                    cf0 = ior(cfa, cfb) 
                    i0 = search_conf(no, nor, lid_0, rid_0, cf0)
                    if (i0 <= 0 .or. i0 > ncf_0) cycle
                    if (cf0 /= conf_0(i0)) cycle
                    g0 = cfgr_0(i0)
                    if (g0 == -1) cycle
                    cfp = 0
                    do o = 0, no - 1
                        if (ibits(cfa, o, 1_8) == 0_8) cycle
                        cfp = ieor(cfp, ibits(cfb, 0_8, o))
                    end do
                    ph = 1
                    do o = 0, no - 1
                        if (ibits(cfp, o, 1_8) == 0_8) ph = -ph
                    end do
                    st_0(g0) = st_0(g0) + st_0(g0) * ph * cffac_0(i0) * cffac_a(ia) * cffac_b(ib)
                end do
            end do
        end do
    end do
end subroutine

end module