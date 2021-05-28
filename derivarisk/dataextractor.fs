namespace DerivaRisk.data

open InstrumentData
open Deedle
open System.Numerics
open MathNet.Numerics.Integration
open MathNet.Numerics.Interpolation
open MathNet.Numerics.LinearAlgebra
open DerivaRisk.Commons

module Extractors=
    type RateIntegrationResults={r_minus_q:float;r:float;q:float}

    let inline private extractor(collection:MDseries array) (id:Identifier):Series<float,float>=
        let spotdata=
            collection
                |> Array.filter(fun x-> x.id = id.id)
                |>Array.toList
               
        match spotdata with
        |[h] -> h.points |> Array.map(fun p-> p.x,p.y) |> Series.ofObservations |> Series.sortByKey
        |_   -> failwith (sprintf "A unique Spot Information for %s is required" (id.id))



    let extract_spot(md:Marketdata)(id:Identifier):Series<float,float>=       
        extractor md.spots id

    let extract_spot_for_time_t(md:Marketdata)(t:float)(id:Identifier):float=       
        let Sx = extract_spot md id
        match (Sx.TryGet(t,Deedle.Lookup.ExactOrGreater),Sx.TryGet(t,Deedle.Lookup.ExactOrSmaller)) with
        |(a,_) when a.HasValue -> a.Value
        |(_,b) when b.HasValue -> b.Value
        |_ -> failwith "No spot data could be retrieved from the provided market data object"                           


    let extract_volatility_surface(md:Marketdata)(id:string)=
         let volsurf = match md.volalities
                            |> Array.filter(fun s -> s.id=id)
                            |> List.ofArray with
                       |[] -> failwith (sprintf "the provided id={%s}" (id))
                       |l -> l
                       

         let spline (smile:Point array) =
            let x = smile |> Array.map(fun p-> p.x)
            let y = smile |> Array.map(fun p-> p.y)
            CubicSpline.InterpolateAkima(x, y);

         let T2smile = volsurf |> List.map(fun v -> v.Maturity,spline v.strike2vol)
                               |> List.sortBy(fun p -> p |> fst)
                               |> Series.ofObservations
         let index = T2smile.Index.KeySequence |> Array.ofSeq
         let volpoint (T:double)(K:double)=
            let t1,t2 = Commons.binSearch T index
            let T1,T2 = index.[t1],index.[t2]
            let v1 = T2smile.Get(T,Lookup.ExactOrSmaller).Interpolate(K)
            let v2 = T2smile.Get(T,Lookup.ExactOrGreater).Interpolate(K)
            if t1=t2 then
                v1
            else
                v1 + (v2-v1)*(T2-T)/(T2-T1)

         volpoint

    let extract_dividends(md:Marketdata)(id:Identifier):Series<float,float>=
        extractor md.dividends id

    let extract_rates(md:Marketdata)(id:Identifier):Series<float,float>=
        extractor md.rates id

    let extract_heston(md:Marketdata)(id:Identifier)=
        let xdata=
            md.stochasticVolatilities
                |> Array.filter(fun x-> x.id = id.id)
                |> Array.toList
               
        match xdata with
        |[h] -> h
        |_   -> failwith (sprintf "A unique Spot Information for %s is required" (id.id))

    let Integrate_Rate(md:Marketdata)(id:Identifier)(t:float,T:float):float=
        let rates_series = extract_rates md id                            
                            |> Series.observations
                            |> Array.ofSeq
                
                    
        let integrand = if rates_series.Length <5 then
                            LinearSpline.Interpolate(rates_series |> Array.map(fst),rates_series |> Array.map(snd)).Interpolate
                         else
                            CubicSpline.InterpolateAkima(rates_series |> Array.map(fst),rates_series |> Array.map(snd)).Interpolate
        
        let rate = GaussLegendreRule.Integrate(System.Func<float,float>(integrand), t, T, 32);
        
        rate/(T-t)


    let Integrate_Dividends(md:Marketdata)(id:Identifier)(t:float,T:float):float=
        
        let dividends_series = extract_dividends md id                            
                                |> Series.observations
                                |> Array.ofSeq

        
        if dividends_series.Length = 0 then
            0.0

        else
            let integrand =
                if dividends_series.Length <5 then
                    LinearSpline.Interpolate(dividends_series |> Array.map(fst),dividends_series |> Array.map(snd)).Interpolate
                else
                    CubicSpline.InterpolateAkima(dividends_series |> Array.map(fst),dividends_series |> Array.map(snd)).Interpolate
            

            
            let div = GaussLegendreRule.Integrate(System.Func<float,float>(integrand), t, T, 32);
            div/(T-t)


    let messageMatrix2Matrix(m:CorrelationMatrix)=
        let v = m.data
        let names = m.assets
        names,
        match m.isRowMajor with
        |true -> Matrix<double>.Build.DenseOfRowMajor(names.Length,names.Length,v)
        | _   -> Matrix<double>.Build.DenseOfColumnMajor(names.Length,names.Length,v)