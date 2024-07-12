import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
import hashlib          # for SHAKE128
import secrets          # for RNG
import numpy as np

# BFV params
deg = 2048
mod = 2**54+1 # CT modulus
mod_t = 5 # Plaintext modulus (prime)

R = polyring_t(deg, mod)
Rt = polyring_t(deg, mod_t) # plaintext

# s = gen_binary_poly(size)
# to do 
s = poly_t.urandom_static(R, 2, secrets.token_bytes(32), 0)
print("BFV.key-gen.sk generated.")
# s.print()
a = poly_t.urandom_static(R, mod, secrets.token_bytes(32),0)
print("BFV.key-gen.a generated.")
# a.print()

# todo: std must be 3.19 (from Illia's repository)
e = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
print("BFV.key-gen.e generated.")
# e.print()

b = -a*s + e
b.print()
pk = (b, a)
print("BFV.key-gen.(pk,sk) generated.")

delta = mod//mod_t
print("BFV.enc.delta generated.")
print(delta)

u = poly_t.urandom_static(R, 2, secrets.token_bytes(32), 0)
print("BFV.enc.u generated.")
# u.print()

e0 = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
e1 = poly_t.grandom_static(R,1, secrets.token_bytes(32), 0)
print("BFV.enc.(e0,e1) generated.")

# -- encrypt
m = poly_t.urandom_static(Rt, 5, secrets.token_bytes(32), 0)
print("BFV.enc.m prepared randomly.")
# m.print()

# ct0 = poly_t(R, None)
# ct1 = poly_t(R, None)
# ct0.set_coeffs()
# ct0 = (pk[0]*u + e0 + delta*m) 
# ct1 = (pk[1]*u + e1). 

# print("m_delta=",m_delta.to_list()[:10])
# print("m*delta=",m.__mul__(delta).to_list()[:10])

m_delta = poly_t(R, (m.__mul__(delta).to_list()))
ct0 = (pk[0]*u + e0 + m_delta)
ct1 = (pk[1]*u + e1)

print("BFV.ct0 =",ct0.to_list()[:10])
print("BFV.ct1 =",ct1.to_list()[:10])

# --- decrypt
s_q = poly_t(R, (s.to_list()))
ct_s = (ct0 + ct1*s_q)
ct_s_coeffs = ct_s.to_list()

ct_s_coeffs = [coeff*mod_t for coeff in ct_s_coeffs]
print("ct_s_coeffs1=",ct_s_coeffs[:10])
ct_s_coeffs = [round(coeff/(mod/2)) for coeff in ct_s_coeffs]
print("ct_s_coeffs2=",ct_s_coeffs[:10])
m_decrypted = [coeff % mod_t for coeff in ct_s_coeffs]

print("m=",m.to_list()[:10])
print("m_decrypted=",m_decrypted[:10])


# b = polyadd(polymul_ntt(negate(a), s, modulus, primitive_root, poly_mod), negate(e), modulus, poly_mod)
# return (b, a), s