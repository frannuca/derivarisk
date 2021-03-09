namespace DerivaRisk

open System
open MathNet.Numerics
open System.Numerics
open MathNet.Numerics.Random
open Deedle

module GBM=
    open MathNet.Numerics.LinearAlgebra
                
     type AssetName = AssetName of string
     let norm = Distributions.Normal()
     let unif = Distributions.ContinuousUniform(0.0,1.0)
     let seedS = 42
     let seedV = 23

     type GBMParams={rate:float;dt:float;sigma:Spot -> float<years> -> float;S0:float}

     let computeMCPaths(mcdata:SimCubeData)(x:GBMParams[])=

         let cube = MCSimCube.generate_cube mcdata
                  
         let rec compute_asset_path_exp(x:GBMParams,S:float[],dZ_asset:float[])(n:int):float[]=
             if(n>=S.Length) then
                 S
             else
                 let T =  float(mcdata.ntimesteps - n)*x.dt * 1.0<years>
                 let sigma = x.sigma (Spot S.[n-1])(T)
                 S.[n]<-S.[n-1]+((x.rate-(sigma**2.0)/2.0)*x.dt+Math.Sqrt(x.dt)*sigma*dZ_asset.[n])
                 compute_asset_path_exp(x,S,dZ_asset)(n+1) 

         
         //cube is a matrix of correlated npaths x assets x time
         let mc = cube
                     |> Array.mapi(fun nsim asset_paths ->                                                        
                             asset_paths
                             |> Array.mapi(fun nasset s_path ->
                                let S = Array.zeroCreate mcdata.ntimesteps
                                S.[0]<-x.[nasset].S0
                                compute_asset_path_exp(x.[nasset],S,s_path)(1)
                             )
                 
                         )
         
         mc,cube 
