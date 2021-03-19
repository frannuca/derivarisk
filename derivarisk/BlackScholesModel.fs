namespace DerivaRisk.model
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Distributions
open DerivaRisk.data
open InstrumentData
open Deedle

type BlackScholesModel()=
    let N x = Normal.CDF(0.0,1.0,x)
    let fN x = Normal.PDF(0.0,1.0,x)

    let validate_data(vanillaoption:OptionBase)=
        if vanillaoption = null then
            failwith "Instrument type must be OptionBase to be priced with Black and Scholes model"

        if vanillaoption.basket.components.Length <> 1 then
            failwith "Vanilla options prices with Black Schole Model must hold only a single underlying"

        vanillaoption

    let mutable greek_delta:Option<float>=None
    let mutable greek_gamma:Option<float>=None
    
    interface IModel with

        member self.npv(instr:Instrument)(md:Marketdata):float=
            
            let vanillaoption = instr :?> OptionBase |> validate_data
            
            let underlyingid = vanillaoption.basket.components.[0]
            let Ss= Extractors.extract_spot md underlyingid
            
            
            let t = float(vanillaoption.startdate.Subtract(md.referenceDate).Hours)/24.0/360.0
            let T = float(vanillaoption.expiry.Subtract(md.referenceDate).Hours)/24.0/360.0
            let K=vanillaoption.strike.amount

            let So=Ss.Get(t,Deedle.Lookup.ExactOrSmaller)
            let sigma_series = Extractors.extract_volatility md underlyingid
            let sigma = sigma_series.Get(t,Deedle.Lookup.ExactOrSmaller)
            let d1 S K T sigma r =
                (log(S/K)+(r+sigma**2.0/2.0)*T)/(sigma*sqrt(T))

            
            let r = Extractors.Integrate_Rate md vanillaoption.currency (t,T)
            let q = Extractors.Integrate_Dividends md vanillaoption.currency (t,T)
            let d_1 = d1 So K (T-t) sigma (r-q)
            let d_2 = d_1-sigma*sqrt(T-t)

            greek_gamma <- Some(fN(d_1)/(So*sigma*sqrt(T-t)))

            let npv,deltag=
                    match vanillaoption.optionType with
                    |OptionType.CALL ->                                         
                                        So*N(d_1)-K*N(d_2),N(d_1)
                
                    |OptionType.PUT ->                                       
                                        K*exp(-(r-q)*(T-t))*N(-d_2)-So*N(-d_1),N(d_1)-1.0
                            
                    |_ -> failwith (sprintf "Unknown Option Type")

            greek_delta <- Some(deltag)
            npv

        member self.greekDelta(instr:Instrument)(md:Marketdata):float array=           

            match greek_delta with
            |Some(x)-> [|x|]
            | _ -> let price = (self :> IModel).npv instr md
                   [|greek_delta |> defaultArg <| 0.0|]
            

        member self.greekGamma(instr:Instrument)(md:Marketdata):float array=
          
            match greek_gamma with
            |Some(x)-> [|x|]
            | _ -> let price = (self :> IModel).npv instr md
                   [|greek_gamma |> defaultArg <| 0.0|]
                    
   