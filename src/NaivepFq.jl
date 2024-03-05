module NaivepFq

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
function RisingPochammer(vs::StaticVector{N,T}) where {N<:Number,T<:Number}
  return RisingPochammer(SVector(vs), (@MVector ones(T, N)))
end
function RisingPochammer(vs::StaticVector{N,T}) where {N<:Number,T}
  return RisingPochammer(SVector(vs), @MVector(T[]))
end

function RisingPochammer(vs::NTuple{N,T}) where {N,T<:Number}
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
Base.:*(p::RisingPochammer, x) = x * prod(p.sₙ)
Base.:/(x, p::RisingPochammer{0}) = x
Base.:/(x, p::RisingPochammer{1}) = x / p.sₙ[1]
Base.:/(x, p::RisingPochammer) = x / prod(p.sₙ)


nontupletypes(αs::A, βs::B, z::Z) where {A<:Tuple{},B<:Tuple{},Z} = Z
nontupletypes(αs::A, βs::B, z::Z) where {A,B<:Tuple{},Z} = promote_type(eltype(A),Z)
nontupletypes(αs::A, βs::B, z::Z) where {A<:Tuple{},B,Z} = promote_type(eltype(B),Z)
nontupletypes(αs::A, βs::B, z::Z) where {A,B,Z} = promote_type(eltype(A),eltype(B),Z)

function nloopsskip(αs, βs, z)
  T = float(nontupletypes(αs, βs, z))
  return max(floor(Int, log2(abs((prod(αs) * z) / prod(βs))) * 16), 1)
end

function pFq(αs, βs, z; maxiters::Int = 1000_000,
             rtol=8eps(float(real(nontupletypes(αs, βs, z)))),
             minloops=nloopsskip(αs, βs, z))
             
  αₙ = RisingPochammer(αs)
  βₙ = RisingPochammer(βs)

  T = float(nontupletypes(αs, βs, z))
  abs(z) < eps(real(T)) && return one(T)

  a, b = one(T), zero(T)
  zⁿ = nbang = one(T)
  n = 0
  while n < minloops || n < maxiters && !fastisapprox(a, b, rtol=rtol, atol=0, nans=true)
    b = a
    zⁿ *= z
    αₙ(n)
    βₙ(n)
    nbang *= (n += 1)
    a += (αₙ * inv(βₙ * nbang)) * zⁿ
    @assert isfinite(a) (αₙ, zⁿ, (βₙ * nbang))
  end
  return a
end


end # module NaivepFq
