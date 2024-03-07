module NaivepFq

# https://arxiv.org/pdf/math/0306302.pdf really useful

using StaticArrays

function fastisapprox(A::Number, B::Number; atol, rtol, nans)
  n = abs2(A - B)
  X, Y = abs2(A), abs2(B)
  (n <= max(abs2(atol), abs2(rtol) * max(X, Y))) && return true
  return nans ? (isnan(X) && isnan(Y)) : false
end

struct RisingPochammer{N,T}
  vs::SVector{N,T}
  sₙ::MVector{N,T}
end
function RisingPochammer(vs::StaticVector{N,T}) where {N,T}
  return RisingPochammer(SVector(vs), (@MVector ones(T, N)))
end
function RisingPochammer(vs::StaticVector{0,T}) where {T}
  return RisingPochammer(SVector(vs), @MVector(T[]))
end

function RisingPochammer(vs::NTuple{N,T}) where {N,T}
  return RisingPochammer(SVector(vs), (@MVector ones(T, N)))
end
#RisingPochammer(vs::Tuple{T}) where {T<:Real} = RisingPochammer(NTuple{1,T}(vs))
function RisingPochammer(vs::NTuple{0})
  RisingPochammer((@SVector Bool[]), (@MVector Bool[]))
end
function RisingPochammer(vs::Tuple{<:Union{}})
  RisingPochammer((@SVector Bool[]), (@MVector Bool[]))
end

(p::RisingPochammer{0})(n) = Returns(nothing)
(p::RisingPochammer)(n) = (p.sₙ .*= p.vs .+ n)
(p::RisingPochammer{0})() = 1
(p::RisingPochammer{1})() = (@inbounds p.sₙ[1])
(p::RisingPochammer)() = prod(p.sₙ)

Base.:*(p::RisingPochammer{0}, x) = x
Base.:*(p::RisingPochammer{1}, x) = p.sₙ[1] * x
Base.:*(p::RisingPochammer, x) = foldl(*, p.sₙ; init=x)
Base.:/(x, p::RisingPochammer{0}) = x
Base.:/(x, p::RisingPochammer{1}) = x / p.sₙ[1]
Base.:/(x, p::RisingPochammer) = foldl(/, p.sₙ; init=x)

nontupletypes(αs::A, βs::B, z::Z) where {A<:Tuple{},B<:Tuple{},Z} = Z
nontupletypes(αs::A, βs::B, z::Z) where {A,B<:Tuple{},Z} = promote_type(eltype(A),Z)
nontupletypes(αs::A, βs::B, z::Z) where {A<:Tuple{},B,Z} = promote_type(eltype(B),Z)
nontupletypes(αs::A, βs::B, z::Z) where {A,B,Z} = promote_type(eltype(A),eltype(B),Z)

function nloopsskip(αs, βs, z)
  T = float(nontupletypes(αs, βs, z))
  val = (foldl(*, αs; init=1) * z) / foldl(*, βs; init=1)
  return max(floor(Int, log2(abs(val)) * 4), 1)
end

#function pFq(αs, βs, z; maxiters::Int = 1000_000,
#             rtol=8eps(float(real(nontupletypes(αs, βs, z)))),
#             minloops=nloopsskip(αs, βs, z))
#             
#  αₙ = RisingPochammer(αs)
#  βₙ = RisingPochammer(βs)
#
#  T = float(nontupletypes(αs, βs, z))
#  abs(z) < eps(real(T)) && return one(T)
#
#  a, b = one(T), zero(T)
#  zⁿ = nbang = one(T)
#  n = 0
#  while n < minloops || n < maxiters && !fastisapprox(a, b, rtol=rtol, atol=0, nans=true)
#    b = a
#    zⁿ *= z
#    αₙ(n)
#    βₙ(n)
#    nbang *= (n += 1)
#    a += (αₙ * inv(βₙ * nbang)) * zⁿ
#  end
#  return isfinite(a) ? a : b
#end

#function pFq(αs, βs, z; maxiters::Int = 1000_000,
#             rtol=8eps(float(real(nontupletypes(αs, βs, z)))),
#             levink=32, levinb=1, levinn=0, minloops=4)
#             
#  αₙ = RisingPochammer(αs)
#  βₙ = RisingPochammer(βs)
#
#  T = float(nontupletypes(αs, βs, z))
#  abs(z) < eps(real(T)) && return one(T)
#
#  a, b = one(T), zero(T)
#  numsum, densum = zero(T), zero(T)
#  s = zⁿ = one(T)
#  jbang = binomialkj = one(float(real(T)))
#  sgn = -1
#  j = 0
#  while j < minloops || j < levink && !fastisapprox(a, b, rtol=rtol, atol=0, nans=true)
#    b = a
#    zⁿ *= z
#    αₙ(j)
#    βₙ(j)
#    jbang *= (j += 1)
#    t = (αₙ * inv(βₙ * jbang)) * zⁿ
#    s += t
#    njb = levinn + j + levinb
#    sgn *= -1
#    den = sgn * binomialkj / t * float(real(T))(njb)^(levink - 2)
#    #den = sgn * binomialkj / t * exp(log(njb)*(levink - 2))
#    num = den * s 
#    densum += den
#    numsum += num
#    a = numsum / densum
#    binomialkj *= (levink - j) / (j + 1)
#  end
#  return a
#end

#function pFq(αs, βs, z; maxiters::Int = 1000_000,
#             rtol=6400eps(float(real(nontupletypes(αs, βs, z)))),
#             minloops=nloopsskip(αs, βs, z))
#             
#  αₙ = RisingPochammer(αs)
#  βₙ = RisingPochammer(βs)
#
#  T = float(nontupletypes(αs, βs, z))
#  abs(z) < eps(real(T)) && return one(T)
#
#  a = one(T)
#  zⁿ = one(typeof(z))
#  nbang = one(float(real(T)))
#  x0 = a
#
#  n = 0
#  zⁿ *= z
#  αₙ(n)
#  βₙ(n)
#  nbang *= (n += 1) # n=1
#  t = (αₙ * inv(βₙ * nbang)) * zⁿ
#  a += t
#  x1 = a
#
#  r = -one(T)
#  p = one(T)
#  q = zero(T)
#
#  while n < minloops || !fastisapprox(r, q, rtol=rtol, atol=0, nans=true)
#    p = q
#    zⁿ *= z
#    αₙ(n)
#    βₙ(n)
#    nbang *= (n += 1)
#    t = (αₙ * inv(βₙ * nbang)) * zⁿ
#    a += t
#    x2 = a
#    isfinite(x2) || return q
#    x1 == x2 && return q
#    den = 1 / (x2 - x1) - 1 / (x1 - x0)
#    iszero(den) && return q
#    q = x1 + 1 / den
#    x0 = x1
#    x1 = x2
#    x0 == x1 && return q
#    r = (1 / n - 1 / (n-1)) / (1 / p / n - 1 / q / (n-1))
#  end
#  return r
#end


#function pFq(αs, βs, z; maxiters::Int = 1000_000,
#             rtol=8eps(float(real(nontupletypes(αs, βs, z)))),
#             minloops=nloopsskip(αs, βs, z))
#             
#  αₙ = RisingPochammer(αs)
#  βₙ = RisingPochammer(βs)
#
#  T = float(nontupletypes(αs, βs, z))
#  abs(z) < eps(real(T)) && return one(T)
#
#  a, b, c = one(T), zero(T), zero(T)
#  zⁿ = nbang = one(T)
#  r, s = zero(T), one(T)
#  invn = one(float(real(T)))
#  n = 0
#  while n < minloops || n < maxiters && !fastisapprox(s, r, rtol=rtol, atol=0, nans=true)
#    c = b
#    b = a
#    zⁿ *= z
#    αₙ(n)
#    βₙ(n)
#    nbang *= (n += 1)
#    t = (αₙ * inv(βₙ * nbang)) * zⁿ
#    isfinite(t) || return a
#    a += t
#    s = r
#    invn_1 = invn
#    invn = 1 / n
#    den = invn * a - invn_1 * b
#    iszero(den) && return r
#    r = a * b * (invn - invn_1) / den
#  end
#  return isfinite(r) ? r : s
#end


function pFq(αs, βs, z; maxiters::Int = 1000_000,
             rtol=8eps(float(real(nontupletypes(αs, βs, z)))),
             drummondk=64, minloops=4)
             
  αₙ = RisingPochammer(αs)
  βₙ = RisingPochammer(βs)

  T = float(nontupletypes(αs, βs, z))
  abs(z) < eps(real(T)) && return one(T)

  a, b = one(T), zero(T)
  numsum, densum = zero(T), zero(T)
  s = zⁿ = one(T)
  jbang = binomialkj = one(float(real(T)))
  sgn = -1
  j = 0
  while j < minloops || j < drummondk && !fastisapprox(a, b, rtol=rtol, atol=0, nans=true)
    b = a
    zⁿ *= z
    αₙ(j)
    βₙ(j)
    jbang *= (j += 1)
    t = (αₙ * inv(βₙ * jbang)) * zⁿ
    s += t
    sgn *= -1
    den = sgn * binomialkj / t
    num = den * s 
    densum += den
    numsum += num
    a = numsum / densum
    binomialkj *= (drummondk - j) / (j + 1)
  end
  return a
end


end # module NaivepFq
