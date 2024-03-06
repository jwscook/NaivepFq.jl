using Test
using NaivepFq, StaticArrays, HypergeometricFunctions

const HFpFq = HypergeometricFunctions.pFq

@testset "NaivepFq" begin
  t0 = 0.0
  t1 = 0.0
  for i in 0:0, j in 1:1, T in (Float64, )
    for k in 0:1, _ in 1:100
      a = Tuple(rand(T, i) .* 10^k)
      b = Tuple(rand(T, j) .* 10^k)
      z = rand(T) .* 10^k
      a1 = SVector(a)
      b1 = SVector(b)
      #try
      t1 += @elapsed r1 = HFpFq(a, b, z)
      t0 += @elapsed r0 = NaivepFq.pFq(a1, b1, z; minloops=10)
      @test r0 ≈ r1
      #catch
      #end
    end
  end
  @show t0, t1
end

#@btime pFq(Tuple(()), (2.0,), 3.0)
#@btime HFpFq(Tuple(()), (2.0,), 3.0)
#@btime pFq(Tuple(()), (2.0+im,), 3.0+im)
#@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im)
#@btime pFq(Tuple(()), (2.0+im,), 3.0+im, rtol=sqrt(eps()))
#@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im, rtol=sqrt(eps()))
#@btime pFq(Tuple(()), (2.0+im,), 3.0+im, rtol=eps())
#@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im, rtol=eps())
#@btime pFq((@SVector Float64[]), (@MVector [2.0]), 3.0)
#@btime HFpFq(Tuple(()), (2.0,), 3.0)
##@btime pFq(Tuple(()) , (2.0,), 3.0)
##@btime pFq(Tuple(()) , (2.0,), 3.0+im)
##@btime pFq(Tuple(()) , (2.0+im,), 3.0)
##@btime pFq(Tuple(()) , (2.0+im,), 3.0+im)
#@btime pFq(Tuple(()), (Dual(2.0+im, 1),), 3.0+im)
#@btime HFpFq(Tuple(()), (Dual(2.0+im, 1),), 3.0+im)
#@btime pFq(Tuple(()), (Dual(20.0+im, 1),), 30.0-im)
#@btime HFpFq(Tuple(()), (Dual(20.0+im, 1),), 30.0-im)
#@btime pFq((1.0, 2.0, 3.0), (Dual(20.0+im, 1),), 30.0-im)
#@btime HFpFq((1.0, 2.0, 3.0), (Dual(20.0+im, 1),), 30.0-im)

