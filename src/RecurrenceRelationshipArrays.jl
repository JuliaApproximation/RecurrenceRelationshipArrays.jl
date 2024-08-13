module RecurrenceRelationshipArrays
using LinearAlgebra, InfiniteArrays, LazyArrays
import LazyArrays: AbstractCachedArray, LazyArrayStyle
import Base: size, getindex, broadcasted, copy, view, +, -
export RecurrenceArray

include("recurrence.jl")

end # module RecurrenceRelationshipArrays
