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

     type GBMParams={rate:float;dt:float;sigma:float;S0:float}

     let computeMCPaths(assets_rho:Matrix<float>,                       
                        npaths:int,
                        ntimesteps:int,
                        x:GBMParams[])=

         let cube = MCSimCube.generate_cube assets_rho seedS npaths ntimesteps (assets_rho.ColumnCount)         
                  
         let rec compute_asset_path_exp(x:GBMParams,S:float[],dZ_asset:float[])(n:int):float[]=
             if(n>=S.Length) then
                 S
             else                 
                 S.[n]<-S.[n-1]+((x.rate-(x.sigma**2.0)/2.0)*x.dt+Math.Sqrt(x.dt)*x.sigma*dZ_asset.[n])
                 compute_asset_path_exp(x,S,dZ_asset)(n+1) 

         
         //cube is a matrix of correlated npaths x assets x time
         let mc = cube
                     |> Array.mapi(fun nsim asset_paths ->                                                        
                             asset_paths
                             |> Array.mapi(fun nasset s_path ->
                                let S = Array.zeroCreate ntimesteps
                                S.[0]<-x.[nasset].S0
                                compute_asset_path_exp(x.[nasset],S,s_path)(1)
                             )
                 
                         )
         
         mc,cube 
