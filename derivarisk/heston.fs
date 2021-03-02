namespace DerivaRisk

open System
open MathNet.Numerics
open System.Numerics
open MathNet.Numerics.Random

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

    type AssetName = AssetName of string
    let norm = Distributions.Normal()
    let unif = Distributions.ContinuousUniform(0.0,1.0)
    let computeMCPaths(assets_rho:Matrix<float>,npaths:int,ntimesteps:int,hestonparams:Map<AssetName,HestonParams>)=
        
        
        let cube = MCSimCube.generate_cube assets_rho 42 npaths ntimesteps (hestonparams.Count)
        let cube_vols = MCSimCube.generate_cube (Matrix<float>.Build.Diagonal(hestonparams.Count,hestonparams.Count,1.0)) 23 npaths ntimesteps (hestonparams.Count)

        let asset_names = hestonparams |> Seq.map(fun kv -> kv.Key) |> Array.ofSeq
        let Ks = hestonparams
                            |> Map.map(fun key v ->                                
                                let K0 = -v.kappa*v.rho*v.theta/v.sigma*v.dt
                                let K1 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_1*v.dt-v.rho/v.sigma
                                let K2 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_2*v.dt+v.rho/v.sigma
                                let K3 = (1.0-v.rho**2.0)*v.gamma_1*v.dt
                                let K4 = (1.0-v.rho**2.0)*v.gamma_2*v.dt
                                let A = K2+0.5*K4
                                {K0=K0;K1=K1;K2=K2;K3=K3;K4=K4;A=A}                            
                            ) 

        let rec compute_asset_path(x:HestonParams,K:Kparams,V:float[],S:float[],dZ_asset:float[],dZ_vol:float[])(n:int):float[]=
            if(n>=S.Length) then
                S
            else                
                let dZ= dZ_asset.[n]
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

                S.[n]<-S.[n-1]*Math.Exp((x.rate-x.dividends)*x.dt+K0+K.K1*V.[n-1]+K.K2*V.[n]+Math.Sqrt(K.K3*V.[n-1]+K.K4*V.[n])*dZ)
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
        
        mc    
        