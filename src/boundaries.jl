abstract type AbstractBoundaryCondition end
struct ZeroPadded <: AbstractBoundaryCondition end
struct Periodic <: AbstractBoundaryCondition end