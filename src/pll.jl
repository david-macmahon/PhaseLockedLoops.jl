using DSP
using Plots

function make_testsig(cf=1/40; N=1000, twin=:, noise_level=0)
    testsig = noise_level * (2*rand(ComplexF32, N) .- (1+1im))
    testsig[twin] .+= cispi.(2*cf .* axes(testsig,1)[twin])
end

function pll(sig; reffreq=1, cutoff=0.25, gain=1e-3)
    responsetype = Lowpass(cutoff)
    designmethod = Butterworth(2)
    lpfbiquad = convert(Biquad, digitalfilter(responsetype, designmethod))
    loop_lpf = DF2TFilter(lpfbiquad, zeros(eltype(sig), 2))

    refsig = one(eltype(sig))
    pll_integral = zero(real(eltype(sig)))
    pllints = zeros(eltype(pll_integral), size(sig))
    refout = zeros(sig)
    pdi = zeros(sig)
    pdf = zeros(sig)
    for (t, x) in enumerate(sig)
        # BEGIN PLL block
        pd_instantaneous = x * conj(refsig)
        push!(pdi, pd_instantaneous)
        pd_filtered = filt(loop_lpf, [pd_instantaneous])[1]
        push!(pdf, pd_filtered)
        pll_integral += angle(pd_filtered) * gain
        push!(pllints, pll_integral)
        refsig = cispi(2 * reffreq * (t + pll_integral))
        push!(refout, refsig)
        # END PLL block
    end

    refout, pdi, pdf, pllints
end