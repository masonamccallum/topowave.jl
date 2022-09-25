# to run test use the function retest() or runtests()
module topowaveTest
using ReTest, UnPack, StartUpDG
using topowave

const_initial(x,y) = reshape([ 1.0+0.0*im 1.0+0.0*im],1,2);
linear_in_x_initial(x,y) = reshape([x + 0.0*im x + 0.0*im],1,2);
linear_in_y_initial(x,y) = reshape([y + 0.0*im y + 0.0*im],1,2);
quadratic_in_x_initial(x,y) = reshape([x^2 + 0.0*im x^2 + 0.0*im],1,2);
quadratic_in_y_initial(x,y) = reshape([y^2 + 0.0*im y^2 + 0.0*im],1,2);

# TODO: add more tests

@testset "test_drs1" begin
    u, parameters = topowave.init_wave(const_initial)
    test_u = topowave.drs1u_test(copy(u), u, parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12

    u, parameters = topowave.init_wave(linear_in_x_initial)
    test_u = topowave.drs1u_test(copy(u),u,parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12

    u, parameters = topowave.init_wave(linear_in_y_initial)
    test_u = topowave.drs1u_test(copy(u),u,parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12
    end

@testset "test_dss2" begin
    u, parameters = topowave.init_wave(const_initial)
    test_u = topowave.dss2u_test(copy(u), u, parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12

    u, parameters = topowave.init_wave(linear_in_x_initial)
    test_u = topowave.dss2u_test(copy(u),u,parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12

    u, parameters = topowave.init_wave(linear_in_y_initial)
    test_u = topowave.dss2u_test(copy(u),u,parameters)
    display(test_u)
    @test test_u ≈ zeros(Complex,size(u,1),size(u,2),size(u,3)) atol=1e-12
    end
end