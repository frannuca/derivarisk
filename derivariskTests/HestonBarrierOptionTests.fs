﻿module HestonMultiAssetBarrierOptionTest
open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open NUnit.Framework
open GBM
open System.IO

[<SetUp>]
let Setup () =
    ()
let nsim=5000
let ntime=100
let T=0.25
let hestonParams = {    HestonParams.dividends=0.02;
                        rate=0.03;
                        kappa=6.2;
                        theta=0.06;
                        sigma=0.5;
                        rho= -0.7;
                        lambda=0.0;
                        S0=100.0;
                        V0=0.03;
                        dt=T/float ntime;
                        gamma_1=0.5;
                        gamma_2=0.5;
                        threshold_phi=1.5}

let duration f = 
     let timer = new System.Diagnostics.Stopwatch()
     timer.Start()
     let returnValue = f()
     printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
     returnValue


[<Test>]
let TestPriceUpKnockInCall()=      
    //let rho = Matrix<float>.Build.DenseOfRowArrays([|1.0|])
    
    let K=90.0
    let BL = 110.0
    let hestonparams=Map([Heston.AssetName("A"), hestonParams])

    let heston_call = Heston.heston_analytical( hestonParams.S0,
                                                K,
                                                hestonParams.kappa,
                                                0.0,
                                                T,
                                                hestonParams.rate,
                                                hestonParams.dividends,
                                                hestonParams.rho,
                                                hestonParams.sigma,
                                                hestonParams.theta,
                                                hestonParams.V0,
                                                true)

 
    let assets = [Heston.AssetName "A"]
    let rho = Matrix<float>.Build.DenseOfRowArrays([|[|1.0|] |]) |> Utils.matrix2Frame assets assets
    let getpar =
        [(Heston.AssetName("A"),hestonParams)]|>Map.ofSeq
    let barrier = UpIn BL
    let mcprice= HestonBarrier.barrierprice(nsim,ntime,getpar,rho,rho,K,barrier,CALL)
    
    printfn "Price Call Option=%f MC price =%f " heston_call mcprice
    Assert.Less(mcprice,heston_call)

