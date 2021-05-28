namespace DerivaRisk.Commons
open System
open MathNet.Numerics.LinearAlgebra
open Deedle

[<AutoOpen>]
module Commons=
    [<Measure>] type years
    [<Measure>] type months
    [<Measure>] type days

    let binSearch (target:float) (arr:float array) =

        let rec ibinsearch (l:int)(h:int) =
            match (h-l) with
                | 0 | 1 -> (l,h)
                | _ ->  let m =  (l+h)/2
                        if arr.[m]>=target then
                            ibinsearch l m
                        else
                            ibinsearch m h
        ibinsearch 0 ((arr |> Array.length) - 1)