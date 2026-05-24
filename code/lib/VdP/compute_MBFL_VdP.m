function M_MBFL = compute_MBFL_VdP()
    M_MBFL.etta =   @(q)(q(1) - 0.5*(1 - q(1)^2)*q(2));
    M_MBFL.gamma =  @(q)(1 - q(2)^2);
    M_MBFL.T =      @(q)(q);
end