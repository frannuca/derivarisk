namespace DerivaRisk.model
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Distributions
open DerivaRisk.data
open InstrumentData
open Deedle

type BlackScholesModel()=
    let N x = Normal.CDF(0.0,1.0,x)

    let validate_data(vanillaoption:OptionBase)=
        if vanillaoption = null then
            failwith "Instrument type must be OptionBase to be priced with Black and Scholes model"

        if vanillaoption.basket.components.Length <> 1 then
            failwith "Vanilla options prices with Black Schole Model must hold only a single underlying"

        vanillaoption

    interface IModel with
        member self.compute (flags:ComputationSelection array)(instr:Instrument)(md:Marketdata):ModelResult=
            
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
                                   
            let npv=
                    match vanillaoption.optionType with
                    |OptionType.CALL -> 
                                        let d_1 = d1 So K (T-t) sigma (r-q)
                                        let d_2 = d_1-sigma*sqrt(T-t)
                                        So*N(d_1)-K*N(d_2)
                
                    |OptionType.PUT ->
                                        let d_1 = d1 So K (T-t) sigma (r-q)
                                        let d_2 = d_1-sigma*sqrt(T)
                                        K*exp(-(r-q)*(T-t))*N(-d_2)-So*N(-d_1)
                            
                    |_ -> failwith (sprintf "Unknown Option Type")

            {ModelResult.npv=npv;delta=None;gamma=None;theta=None;rho=None}