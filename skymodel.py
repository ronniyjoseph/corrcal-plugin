def sky_moment(n_order, k1=4100, gamma1=1.59, k2=4100, gamma2=2.5, s_low=1e-5, s_mid=1, s_high=10.):
    # Check whether the breakpoints are correct
    if s_high < s_mid:
        s_mid = s_high
    if s_low > s_mid:
        s_mid = s_low

    moment = k1 / (n_order + 1 - gamma1) * (s_mid ** (n_order + 1 - gamma1) - s_low ** (n_order + 1 - gamma1)) + \
             k2 / (n_order + 1 - gamma2) * (s_high ** (n_order + 1 - gamma2) - s_mid ** (n_order + 1 - gamma2))

    return moment
