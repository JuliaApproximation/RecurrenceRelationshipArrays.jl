module RecurrenceRelationshipArrays
using LinearAlgebra, InfiniteArrays, LazyArrays, FillArrays, ArrayLayouts, BandedMatrices
import LazyArrays: AbstractCachedArray, LazyArrayStyle, AbstractPaddedLayout, resizedata!, zero!
import Base: axes, size, getindex, broadcasted, copy, view, +, -, Slice, tail
import FillArrays: AbstractFill
import BandedMatrices: AbstractBandedMatrix
import InfiniteArrays: OneToInf, AbstractInfUnitRange
import ArrayLayouts: MatMulVecAdd
import RecurrenceRelationships: forwardrecurrence_next, forwardrecurrence_partial!
export RecurrenceArray

const LazyArraysBandedMatricesExt = Base.get_extension(LazyArrays, :LazyArraysBandedMatricesExt)
const AbstractLazyBandedLayout = LazyArraysBandedMatricesExt.AbstractLazyBandedLayout

include("clenshaw.jl")
include("recurrence.jl")

end # module RecurrenceRelationshipArrays
