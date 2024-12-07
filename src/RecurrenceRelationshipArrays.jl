module RecurrenceRelationshipArrays
using LinearAlgebra, InfiniteArrays, LazyArrays, FillArrays, ArrayLayouts, BandedMatrices
import LazyArrays: AbstractCachedArray, LazyArrayStyle, AbstractPaddedLayout, resizedata!, zero!, paddeddata
import Base: axes, size, getindex, broadcasted, copy, view, +, -, Slice, tail
import FillArrays: AbstractFill, getindex_value
import BandedMatrices: AbstractBandedMatrix, bandwidths, _BandedMatrix
import InfiniteArrays: OneToInf, AbstractInfUnitRange
import ArrayLayouts: MatMulVecAdd, sublayout, MemoryLayout, sub_materialize, transposelayout, materialize!, _fill_lmul!
import RecurrenceRelationships: forwardrecurrence_next, forwardrecurrence_partial!, check_clenshaw_recurrences, clenshaw, polynomialtype
export RecurrenceArray, Clenshaw

const LazyArraysBandedMatricesExt = Base.get_extension(LazyArrays, :LazyArraysBandedMatricesExt)
const AbstractLazyBandedLayout = LazyArraysBandedMatricesExt.AbstractLazyBandedLayout
const LazyBandedLayout = LazyArraysBandedMatricesExt.LazyBandedLayout


include("clenshaw.jl")
include("recurrence.jl")

end # module RecurrenceRelationshipArrays
