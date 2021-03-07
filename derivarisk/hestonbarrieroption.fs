namespace DerivaRisk
open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open System.IO
open Deedle
module HestonBarrier=
    
    type BARRIERTYPE=
        |DownOut of float
        |UpOut of float
        |DownIn of float
        |UpIn of float

    let duration f = 
         let timer = new System.Diagnostics.Stopwatch()
         timer.Start()
         let returnValue = f()
         printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
         returnValue


    let barrierprice(nsim:int,
                             ntime:int,
                             hestonParams:Map<Heston.AssetName,Heston.HestonParams>,
                             rhoS:Frame<Heston.AssetName,Heston.AssetName>,
                             rhoVol:Frame<Heston.AssetName,Heston.AssetName>,
                             K:float,
                             barrier:BARRIERTYPE,
                             optionType:Utils.OptionType
                             )=

        let mcpaths, _ = Heston.computeMCPaths(rhoS,rhoVol,nsim,ntime,hestonParams)

        let terminalpayoff S =
            match optionType with
            |Utils.CALL -> Math.Max(S-K,0.0)
            |Utils.PUT-> Math.Max(K-S,0.0)

        let tpayoff= match barrier with
                        |UpIn(bl) ->
                            let payoff path =
                                let internal_payoff active spot=
                                    if spot >= bl || active then
                                        terminalpayoff spot, true
                                    else 
                                        0.0, false
        
                                let (sT,active) = path |> Array.fold(fun active s -> internal_payoff (active|>snd) s) (0.0,false)
                                sT
                            payoff
                        |DownIn(bl) ->
                            let payoff path =
                                let internal_payoff active spot=
                                    if spot <= bl || active then
                                        terminalpayoff spot, true
                                    else 
                                        0.0, false
        
                                let (sT,active) = path |> Array.fold(fun active s -> internal_payoff (active|>snd) s) (0.0,false)
                                sT
                            payoff
                        |UpOut(bl) ->
                            let payoff path =
                                let internal_payoff knocked spot=
                                    if spot >= bl || knocked then
                                        0.0,true
                                    else 
                                        terminalpayoff spot, false
        
                                let (sT,active) = path |> Array.fold(fun active s -> internal_payoff (active|>snd) s) (0.0,false)
                                sT
                            payoff
                        |DownOut(bl) ->
                            let payoff path =
                                let internal_payoff knocked spot=
                                    if spot <= bl || knocked then
                                        0.0,true
                                    else 
                                        terminalpayoff spot, false
        
                                let (sT,active) = path |> Array.fold(fun active s -> internal_payoff (active|>snd) s) (0.0,false)
                                sT
                            payoff
                        

        mcpaths
        |> Array.map(fun assets -> assets|> Array.map(tpayoff) |> Array.sum)
        |> Array.average
