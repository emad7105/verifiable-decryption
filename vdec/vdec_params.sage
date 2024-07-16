vname = "param"               # variable name

deg   = 2048                  # ring Rp degree d
mod   = 2**54+1                 # CT modulus
mod_t = 5 # Plaintext modulus (prime)
dim   = (2,3)                 # dimensions of A in Rp^(m,n)

v_e = 1

# B_e
# sigma_err = 1.55*2 # this set in grandom_static since we used log2o=1
# B_e = sqrt(2*deg)*sigma_err
# print("B_e = ", B_e)

# B_v
# delta_phi_m = deg # to-do check!!?
# inf_norm_t = mod_t
# B_v = mod / (2*delta_phi_m*inf_norm_t) - 0.5
# print("B_v = ", B_v)

#       [ s,    e,  v_inhv]
wpart = [ [0], [1], [2] ]   # partition of s
wl2   = [ 0 , sqrt(2*deg)*(1.55*2), mod / (2*deg*mod_t) - 0.5  ]  # l2-norm bounds: l2(s) <= sqrt(2048)
wbin  = [ 1, 0, 0       ]  # binary coeffs
wrej  = [ 0             ]  # rejection sampling