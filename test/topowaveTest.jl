# to run test use the function retest() or runtests()
module topowaveTest
using ReTest, topowave
x=1
@testset "test group 1" begin
    @test x==1
end

@testset "test group 2" begin
    @test sim()!=0
end
end