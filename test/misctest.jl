struct TestModel <: AbstractModel end

tm = TestModel()

gillespie(tm, 10)
