namespace DerivaRisk

open System
open MathNet.Numerics
open System.Numerics
open MathNet.Numerics.Random
open Deedle
open System.Numerics
open MathNet.Numerics.Integration

//https://quant.stackexchange.com/questions/18684/heston-model-option-price-formula

 module Heston=
    open MathNet.Numerics.LinearAlgebra

    type HestonParams={
        S0:float;
        V0:float;
        kappa:float;
        lambda:float;        
        rate:float;
        dividends:float;
        rho:float;
        sigma:float;
        theta:float;
        dt:float;
        gamma_1:float;
        gamma_2:float;
        threshold_phi:float
    } with member self.E = Math.Exp(-self.kappa*self.dt)

    type Kparams={K0:float;K1:float;K2:float;K3:float;K4:float;A:float}

    let computeK(v:HestonParams)=
        let K0 = -v.kappa*v.rho*v.theta/v.sigma*v.dt
        let K1 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_1*v.dt-v.rho/v.sigma
        let K2 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_2*v.dt+v.rho/v.sigma
        let K3 = (1.0-v.rho**2.0)*v.gamma_1*v.dt
        let K4 = (1.0-v.rho**2.0)*v.gamma_2*v.dt
        let A = K2+0.5*K4
        {K0=K0;K1=K1;K2=K2;K3=K3;K4=K4;A=A}
        
    type AssetName = AssetName of string    
    let seedS = 42
    let seedV = 1


    let computeMCPaths(assets_rho_frame:Frame<AssetName,AssetName>,
                       vol_rho_frame:Frame<AssetName,AssetName>, 
                       npaths:int,
                       ntimesteps:int,
                       hestonparams:Map<AssetName,HestonParams>)=

        let asset_names = hestonparams |> Seq.map(fun kv -> kv.Key) |> Array.ofSeq
        let assets_rho:MathNet.Numerics.LinearAlgebra.Matrix<float> = Utils.frame2Matrix asset_names asset_names assets_rho_frame
        let vol_rho:MathNet.Numerics.LinearAlgebra.Matrix<float>    = Utils.frame2Matrix asset_names asset_names vol_rho_frame
        
        let cube = MCSimCube.generate_cube assets_rho seedS npaths ntimesteps (hestonparams.Count)
        let cube_vols = MCSimCube.generate_cube vol_rho seedV npaths ntimesteps (hestonparams.Count)

        
        let Ks = hestonparams |> Map.map(fun key v -> computeK v) 
        
        let rec compute_asset_path(x:HestonParams,K:Kparams,V:float[],S:float[],dZ_asset:float[],dZ_vol:float[])(n:int):float[]=
            if(n>=S.Length) then
                S
            else                
                let dZs= dZ_asset.[n]
                let dZv = dZ_vol.[n]
                let m=x.theta + (V.[n-1]-x.theta)*x.E
                let m2 = m**2.0
                let s2 = V.[n-1]*(x.sigma**2.0)*x.E/x.kappa*(1.0-x.E)+ x.theta*(x.sigma**2.0)/(2.0*x.kappa)*(1.0-x.E)**2.0
                let phi = s2/m2
                let K0 =
                        if phi <= x.threshold_phi then
                            let b= Math.Sqrt(2.0/phi-1.0+Math.Sqrt(2.0/phi*(2.0/phi-1.0)))
                            let a= m/(1.0+b**2.0)
                            
                            V.[n] <- a*(b+dZv)**2.0
                            if K.A<1.0/(2.0*a) then
                                let mx=Math.Exp(K.A*b**2.0*a/(1.0-2.0*K.A*a))/Math.Sqrt(1.0-2.0*K.A*a)
                                -Math.Log(mx)-(K.K1+0.5*K.K3)*V.[n-1]
                            else
                                K.K0
                        else
                            let Uv = Distributions.Normal.CDF(0.0,1.0,dZv)
                            let p = (phi - 1.0)/(phi+1.0)
                            let beta = (1.0-p)/m
                            V.[n] <-
                                    if Uv <= p then
                                        0.0
                                    else
                                        1.0/beta*Math.Log((1.0-p)/(1.0-Uv))

                            if K.A<beta then
                                let mx=p+beta*(1.0-p)/(beta-K.A)
                                -Math.Log(mx)-(K.K1+0.5*K.K3)*V.[n-1]
                            else
                                K.K0
                
                S.[n]<-S.[n-1]*Math.Exp((x.rate-x.dividends)*x.dt+K0+K.K1*V.[n-1]+K.K2*V.[n]+Math.Sqrt(K.K3*V.[n-1]+K.K4*V.[n])*dZs)
                compute_asset_path(x,K,V,S,dZ_asset,dZ_vol)(n+1)                                       
        //cube is a matrix of correlated npaths x assets x time
        let mc = cube
                    |> Array.mapi(fun nsim asset_paths ->
                                                       
                            asset_paths
                            |> Array.mapi(fun nasset s_path ->
                                let asset = asset_names.[nasset]
                                let heston_par = hestonparams.[asset]
                                let Kparam = Ks.[asset]
                                let V = Array.zeroCreate(ntimesteps)
                                let S = Array.zeroCreate(ntimesteps)
                                V.[0]<-heston_par.V0
                                S.[0]<-heston_par.S0                                
                                compute_asset_path(heston_par,Kparam,V,S,s_path,cube_vols.[nsim].[nasset])(1)
                                
                            )
                
                        )
        
        mc, cube    


    let heston_analytical(S:float,K:float,kappa:float,l:float,tau:float,r:float,q:float,rho:float,sigma:float,theta:float,vo:float,call:bool)=

        let integrand(phi)=
            let b1=Complex(kappa+l-rho*sigma,0.0)
            let b2 = Complex(kappa + l,0.0)
            let u1=Complex(0.5,0.0)
            let u2= Complex(-0.5,0.0)
            let a = Complex(kappa*theta,0.0)
            let i= Complex(0.0,1.0)
            let sigma2 = Complex(sigma**2.0,0.0)
            let logK = Math.Log(K)
            let logS = Math.Log(S)
            let tau = Complex(tau,0.0)

            let z:Complex = rho*sigma*phi*i

            let g(b:Complex,d:Complex)=
                (b-z-d)/(b-z+d)

            let d(b:Complex,u:Complex)=
                Complex.Sqrt(Complex.Pow(b-z,2.0) - sigma2*(2.0*i*u*phi-phi**2.0))

            let C(b:Complex,d:Complex,g:Complex):Complex=
                let ratio = (1.0-g*Complex.Exp(-d*tau))/(1.0 - g)
                (r-q)*i*phi*tau + a/sigma2 * ((b-z-d)*tau-2.0*Complex.Log(ratio))

            let D(b:Complex,d:Complex,g:Complex)=
                (b-z-d)/(sigma2)*(1.0-Complex.Exp(-d*tau))/(1.0-g*Complex.Exp(-d*tau))

            let d1 = d(b1,u1)
            let d2 = d(b2,u2)
            let g1 = g(b1,d1)
            let g2 = g(b2,d2)
            let C1 = C(b1,d1,g1)
            let C2 = C(b2,d2,g2)
            let D1 = D(b1,d1,g1)
            let D2 = D(b2,d2,g2)
            
            let f1 = (Complex.Exp(C1+D1*vo+i*phi*logS))
            let f2 = (Complex.Exp(C2+D2*vo+i*phi*logS))
            
            let Ia = Complex.Exp(-i*phi*logK)/(i*phi)*(S*Complex.Exp(-q*tau)*f1-K*Complex.Exp(-r*tau)*f2)

            Ia.Real


        let I= GaussLegendreRule.Integrate(Func<float,float>(integrand), 0.0, 100.0, 32);
        let c = 0.5*S*Math.Exp(-q*tau)-0.5*K*Math.Exp(-r*tau)+1.0/Math.PI*I
        if call then                       
            c
        else
            c+K*Math.Exp(-r*tau)-S*Math.Exp(-q*tau)
    



