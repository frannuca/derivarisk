namespace DerivaRisk

module Black=
    open MathNet.Numerics.LinearAlgebra
    open MathNet.Numerics.Distributions

    type BlackData={T:float;S:float;K:float;r:float;sigma:float}

    let N x = Normal.CDF(0.0,1.0,x)

    let d1 S K T sigma r =
        (log(S/K)+(r+sigma**2.0/2.0)*T)/(sigma*sqrt(T))

    let call(x:BlackData):float =
        let d_1 = d1 x.S x.K x.T x.sigma x.r        
        let d_2 = d_1-x.sigma*sqrt(x.T)
        x.S*N(d_1)-x.K*N(d_2)

    let put(x:BlackData):float =
        let d_1 = d1 x.S x.K x.T x.sigma x.r        
        let d_2 = d_1-x.sigma*sqrt(x.T)
        x.K*exp(-x.r*x.T)*N(-d_2)-x.S*N(-d_1)

