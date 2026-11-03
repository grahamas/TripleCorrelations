#= /home/graham/git/TripleCorrelations/src/sequence_class.jl:19 =# LoopVectorization.@turbo_debug(for p0 = eachindex(data)
        #= /home/graham/git/TripleCorrelations/src/sequence_class.jl:20 =#
        contribution += data[p0] * data_λ₁[p0] * data_λ₂[p0]
    end) = begin
    begin
        (var"#loop#bounds#", var"#vargs#") = var"#lv#tuple#args#"
        var"##p0_loop#4536" = Core.getfield(var"#loop#bounds#", 1, false)
        var"##p0_looplen#4537" = ArrayInterface.static_length(var"##p0_loop#4536")
        var"##p0_loopstop#4538" = (rem)(LoopVectorization.last(var"##p0_loop#4536"), Int64)
        VectorizationBase.assume(var"##p0_loopstop#4538" ≥ static(0))
        (var"##vptr##_data", var"##vptr##_data_λ₁", var"##vptr##_data_λ₂") = VectorizationBase.stridedpointers(getfield(var"#vargs#", 1, false))
        var"####op#4543_#or#init#type#" = Core.getfield(var"#vargs#", 2, false)
        var"##Tloopeltype##" = LoopVectorization.promote_type(LoopVectorization.eltype(var"##vptr##_data"), LoopVectorization.eltype(var"##vptr##_data_λ₁"), LoopVectorization.eltype(var"##vptr##_data_λ₂"), var"####op#4543_#or#init#type#")
        var"##Wvecwidth##" = LoopVectorization.pick_vector_width(var"##Tloopeltype##")
        local var"####op#4543__2"
    end
    begin
        let var"##vptr##_data" = LoopVectorization.gesp(var"##vptr##_data", (Static.Zero(),)), var"##vptr##_data#offset##" = VectorizationBase.increment_ptr(var"##vptr##_data"), var"##vptr##_data_λ₁" = LoopVectorization.gesp(var"##vptr##_data_λ₁", (Static.Zero(),)), var"##vptr##_data_λ₁#offset##" = VectorizationBase.increment_ptr(var"##vptr##_data_λ₁"), var"##vptr##_data_λ₂" = LoopVectorization.gesp(var"##vptr##_data_λ₂", (Static.Zero(),)), var"##vptr##_data_λ₂#offset##" = VectorizationBase.increment_ptr(var"##vptr##_data_λ₂")
            let p0 = 0
                (var"####op#4543__2",) = if var"##p0_loopstop#4538" ≥ LoopVectorization.vsub_nsw(LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(4)), static(1))
                        let var"####op#4543__4" = LoopVectorization.zero_vecunroll(static(4), var"##Wvecwidth##", LoopVectorization.of_same_size(var"####op#4543_#or#init#type#", var"##Tloopeltype##"), static(32))
                            begin
                                var"##maybeterm#4545" = true
                                while var"##maybeterm#4545"
                                    var"####op#4539__4" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data", var"##vptr##_data#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 4, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                                    var"####op#4540__4" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₁", var"##vptr##_data_λ₁#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 4, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                                    var"####op#4541__4" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₂", var"##vptr##_data_λ₂#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 4, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                                    var"####op#4542__4" = LoopVectorization.mul_fast(var"####op#4540__4", var"####op#4541__4")
                                    var"####op#4543__4" = LoopVectorization.vfmadd_fast(var"####op#4539__4", var"####op#4542__4", var"####op#4543__4")
                                    p0 = LoopVectorization.vadd_nsw(LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(4)), p0)
                                    var"##maybeterm#4545" = LoopVectorization.vcmpendzs(p0, var"##p0_loop#4536", LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(4)))
                                end
                            end
                            (LoopVectorization.contract_add(var"####op#4543__4", static(2)),)
                        end
                    else
                        var"####op#4543__2" = LoopVectorization.zero_vecunroll(static(2), var"##Wvecwidth##", LoopVectorization.of_same_size(var"####op#4543_#or#init#type#", var"##Tloopeltype##"), static(32))
                        (var"####op#4543__2",)
                    end
                if LoopVectorization.vcmpendzs(p0, var"##p0_loop#4536", LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(2)))
                    var"####op#4539__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data", var"##vptr##_data#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                    var"####op#4540__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₁", var"##vptr##_data_λ₁#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                    var"####op#4541__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₂", var"##vptr##_data_λ₂#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000000, 1}((p0,)), LoopVectorization.False(), static(32))
                    var"####op#4542__2" = LoopVectorization.mul_fast(var"####op#4540__2", var"####op#4541__2")
                    var"####op#4543__2" = LoopVectorization.partialmap(LoopVectorization.vfmadd_fast, var"####op#4543__2", static(2), static(0), var"####op#4539__2", var"####op#4542__2", var"####op#4543__2")
                    p0 = LoopVectorization.vadd_nsw(LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(2)), p0)
                end
                if LoopVectorization.cmpend(p0, var"##p0_loop#4536")
                    let var"##mask##" = LoopVectorization.mask(var"##Wvecwidth##", var"##p0_looplen#4537")
                        if p0 > LoopVectorization.vsub_nsw(var"##p0_loopstop#4538", LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(1)))
                            var"####op#4539__1" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data", var"##vptr##_data#offset##"), (LoopVectorization.MM(var"##Wvecwidth##", p0, static(1)),), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4540__1" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₁", var"##vptr##_data_λ₁#offset##"), (LoopVectorization.MM(var"##Wvecwidth##", p0, static(1)),), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4541__1" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₂", var"##vptr##_data_λ₂#offset##"), (LoopVectorization.MM(var"##Wvecwidth##", p0, static(1)),), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4542__1" = LoopVectorization.mul_fast(var"####op#4540__1", var"####op#4541__1")
                            var"####op#4543__2" = LoopVectorization.ifelselast(LoopVectorization.vfmadd_fast, var"##mask##", static(1), static(0), var"####op#4539__1", var"####op#4542__1", var"####op#4543__2")
                            p0 = LoopVectorization.vadd_nsw(var"##Wvecwidth##", p0)
                        else
                            var"####op#4539__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data", var"##vptr##_data#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000002, 1}((p0,)), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4540__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₁", var"##vptr##_data_λ₁#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000002, 1}((p0,)), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4541__2" = LoopVectorization._vload(VectorizationBase.reconstruct_ptr(var"##vptr##_data_λ₂", var"##vptr##_data_λ₂#offset##"), LoopVectorization.Unroll{1, Int(var"##Wvecwidth##"), 2, 1, Int(var"##Wvecwidth##"), 0x0000000000000002, 1}((p0,)), var"##mask##", LoopVectorization.False(), static(32))
                            var"####op#4542__2" = LoopVectorization.mul_fast(var"####op#4540__2", var"####op#4541__2")
                            var"####op#4543__2" = LoopVectorization.ifelselast(LoopVectorization.vfmadd_fast, var"##mask##", static(2), static(0), var"####op#4539__2", var"####op#4542__2", var"####op#4543__2")
                            p0 = LoopVectorization.vadd_nsw(LoopVectorization.vmul_nsw(var"##Wvecwidth##", static(2)), p0)
                        end
                    end
                end
            end
        end
        var"####op#4543_##onevec##" = LoopVectorization.collapse_add(var"####op#4543__2")
    end
end