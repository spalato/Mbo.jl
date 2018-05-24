export linear, R1, R2, R3, R4

#=                  RESPONSES                                           =#
#=                  first order                                         =#
linear(t, ω_eg, g) = cis(-ω_eg*t) * exp(-g(t))

#=                  third order                                         =#
# TODO: check if defining conjugates in function is more costly than passing
# in. Should not due to dispatch
@inline function F1_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    g_cc′(t) = (conj ∘ g_cc)(t)
    g_bb′(t) = (conj ∘ g_bb)(t)
    g_ba′(t) = (conj ∘ g_ba)(t)
    g_ca′(t) = (conj ∘ g_ca)(t)
    g_cb′(t) = (conj ∘ g_cb)(t)
    (
        - g_cc′(t2)
        - g_bb′(t3)
        - g_aa(t1+t2+t3)
        - g_cb′(t2+t3)
        + g_cb′(t2)
        + g_cb′(t3)
        + g_ca(t1+t2)
        - g_ca(t1)
        + g_ca′(t2+t3)
        - g_ca′(t3)
        + g_ba(t1+t2+t3)
        - g_ba(t1+t2)
        + g_ba′(t3)
    )
end
@inline function F2_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    g_cc′(t) = (conj ∘ g_cc)(t)
    g_bb′(t) = (conj ∘ g_bb)(t)
    g_ba′(t) = (conj ∘ g_ba)(t)
    g_ca′(t) = (conj ∘ g_ca)(t)
    g_cb′(t) = (conj ∘ g_cb)(t)
    (
        - g_cc′(t1+t2)
        - g_bb′(t3)
        - g_aa(t2+t3)
        - g_cb′(t1+t2+t3)
        + g_cb′(t1+t2)
        + g_cb′(t3)
        + g_ca(t2)
        + g_ca′(t1+t2+t3)
        - g_ca′(t1)
        - g_ca′(t3)
        + g_ba(t2+t3)
        - g_ba(t2)
        + g_ba′(t3)
    )
end
@inline function F3_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    g_cc′(t) = (conj ∘ g_cc)(t)
    g_bb′(t) = (conj ∘ g_bb)(t)
    g_ba′(t) = (conj ∘ g_ba)(t)
    g_ca′(t) = (conj ∘ g_ca)(t)
    g_cb′(t) = (conj ∘ g_cb)(t)
    (
        - g_cc′(t1)
        - g_bb′(t2+t3)
        - g_aa(t3)
        - g_cb′(t1+t2+t3)
        + g_cb′(t1)
        + g_cb′(t2+t3)
        + g_ca′(t1+t2+t3)
        - g_ca′(t1+t2)
        - g_ca′(t2+t3)
        + g_ca′(t2)
        + g_ba(t3)
        - g_ba′(t2+t3)
        + g_ba′(t2)
    )
end
@inline function F4_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    g_cc′(t) = (conj ∘ g_cc)(t)
    g_bb′(t) = (conj ∘ g_bb)(t)
    g_ba′(t) = (conj ∘ g_ba)(t)
    g_ca′(t) = (conj ∘ g_ca)(t)
    g_cb′(t) = (conj ∘ g_cb)(t)
    (
        - g_cc(t3)
        - g_bb(t2)
        - g_aa(t1)
        - g_cb(t2+t3)
        + g_cb(t2)
        + g_cb(t3)
        - g_ca(t1+t2+t3)
        + g_ca(t1+t2)
        + g_ca(t2+t3)
        - g_ca(t2)
        - g_ba(t1+t2)
        + g_ba(t1)
        + g_ba(t2)
    )
end

#= Correponds to the 4 point correlation function: <V_ga V_ab V_bc V_cg>
    and the Feynman Diagram:
      --|c g|
      --|b g|
      --|a g|
 =#    
"""
    R1(t1, t2, t3, wa, wb, wc, lineshapes...)

Compute R1 response for the Hilbert cycle <g a b c g>.
"""
@inline function R1(t1, t2, t3, wa, wb, wc, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    wab = wa-wb
    wac = wa-wc
    (cis(-wab*t3-wac*t2-wa*t1)
        *exp(F1_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)))
end
@inline function R2(t1, t2, t3, wa, wb, wc, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    wab = wa-wb
    wac = wa-wc
    (cis(-wab*t3-wac*t2+wc*t1)
        *exp(F2_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)))
end
@inline function R3(t1, t2, t3, wa, wb, wc, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    wab = wa-wb
    wac = wa-wc
    (cis(-wab*t3+wb*t2+wc*t1)
        *exp(F3_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)))
end
@inline function R4(t1, t2, t3, wa, wb, wc, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)
    (cis(-wc*t3-wb*t2-wa*t1)
        *exp(F4_arg(t1, t2, t3, g_aa, g_bb, g_cc, g_ba, g_ca, g_cb)))
end
