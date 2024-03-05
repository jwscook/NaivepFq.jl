
using NaivepFq, StaticArrays, HypergeometricFunctions

const HFpFq = HypergeometricFunctions.pFq

@btime pFq(Tuple(()), (2.0,), 3.0)
@btime HFpFq(Tuple(()), (2.0,), 3.0)
@btime pFq(Tuple(()), (2.0+im,), 3.0+im)
@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im)
@btime pFq(Tuple(()), (2.0+im,), 3.0+im, rtol=sqrt(eps()))
@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im, rtol=sqrt(eps()))
@btime pFq(Tuple(()), (2.0+im,), 3.0+im, rtol=eps())
@btime HFpFq(Tuple(()), (2.0+im,), 3.0 + im, rtol=eps())
@btime pFq((@SVector Float64[]), (@MVector [2.0]), 3.0)
@btime HFpFq(Tuple(()), (2.0,), 3.0)
#@btime pFq(Tuple(()) , (2.0,), 3.0)
#@btime pFq(Tuple(()) , (2.0,), 3.0+im)
#@btime pFq(Tuple(()) , (2.0+im,), 3.0)
#@btime pFq(Tuple(()) , (2.0+im,), 3.0+im)
@btime pFq(Tuple(()), (Dual(2.0+im, 1),), 3.0+im)
@btime HFpFq(Tuple(()), (Dual(2.0+im, 1),), 3.0+im)
@btime pFq(Tuple(()), (Dual(20.0+im, 1),), 30.0-im)
@btime HFpFq(Tuple(()), (Dual(20.0+im, 1),), 30.0-im)
@btime pFq((1.0, 2.0, 3.0), (Dual(20.0+im, 1),), 30.0-im)
@btime HFpFq((1.0, 2.0, 3.0), (Dual(20.0+im, 1),), 30.0-im)
u
