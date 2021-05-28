namespace DerivaRisk.model
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Distributions
open DerivaRisk.data
open InstrumentData
open Deedle
open System
open MathNet.Numerics
open System.Numerics
open MathNet.Numerics.Random
open Deedle
open MathNet.Numerics.Integration
open InstrumentData

type HestonAnalytical()=
    
    let validate_data(vanillaoption:OptionBase)=
        if vanillaoption = null then
            failwith "Instrument type must be OptionBase to be priced with Black and Scholes model"

        if vanillaoption.basket.components.Length <> 1 then
            failwith "Vanilla options prices with Black Schole Model must hold only a single underlying"

        vanillaoption


    interface IModel with
        member self.npv(instr:Instrument)
                       (md:Marketdata):float=

            let vanillaoption = instr :?> OptionBase |> validate_data

            let underlyingId = vanillaoption.basket.components.[0]

            let t = float(vanillaoption.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            let T = float(vanillaoption.expiry.Subtract(md.referenceDate).TotalDays)/360.0

            let S = Extractors.extract_spot_for_time_t md t underlyingId 

            let K= vanillaoption.strike.amount

            let hestonparam = Extractors.extract_heston md underlyingId

            let kappa = hestonparam.kappa
            let l = hestonparam.lambda
            let tau = T
            let rho = hestonparam.rho
            let sigma = hestonparam.sigma
            let theta = hestonparam.theta
            let vo = hestonparam.Vo
           
                                   
            let r = Extractors.Integrate_Rate md vanillaoption.currency (t,T)
            let q = Extractors.Integrate_Dividends md vanillaoption.currency (t,T)
            let r_minus_q = r-q

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
                    (r_minus_q)*i*phi*tau + a/sigma2 * ((b-z-d)*tau-2.0*Complex.Log(ratio))

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
            let npv=
                    match vanillaoption.optionType with                       
                    | OptionType.CALL -> c
                    | _ -> c+K*Math.Exp(-r*tau)-S*Math.Exp(-q*tau)
            npv

        member self.greekDelta(instr:Instrument)
                              (md:Marketdata):float array=

            let model = self:>IModel
            let vanillaoption = instr:?>OptionBase
            let npv = model.npv instr md
            let md_p = deepcopy.shockMarketData(md,0.01)
            let npv_p = model.npv instr md_p

            let t = float(vanillaoption.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            let T = float(vanillaoption.expiry.Subtract(md.referenceDate).TotalDays)/360.0

            let vanillaoption = instr :?> OptionBase |> validate_data            
            let underlyingId = vanillaoption.basket.components.[0]
            let So = Extractors.extract_spot_for_time_t md t underlyingId 
            let Sp = Extractors.extract_spot_for_time_t md_p t underlyingId
            let deltagreek = (npv_p-npv)/(Sp-So)
            [|deltagreek|]


        member self.greekGamma(instr:Instrument)
                              (md:Marketdata):float array=

            let model = self:>IModel
            let vanillaoption = instr:?>OptionBase
            let npv = model.npv instr md
            let h=0.01

            let md_hp = deepcopy.shockMarketData(md,h)
            let npv_hp = model.npv instr md_hp

            let md_hm = deepcopy.shockMarketData(md,-h)
            let npv_hm = model.npv instr md_hm


            let t = float(vanillaoption.startdate.Subtract(md.referenceDate).TotalDays)/360.0            

            let vanillaoption = instr :?> OptionBase |> validate_data            
            let underlyingId = vanillaoption.basket.components.[0]
            let So = Extractors.extract_spot_for_time_t md t underlyingId 
            let Sp = Extractors.extract_spot_for_time_t md_hp t underlyingId

            let deltagreek = (npv_hp-2.0*npv+npv_hm)/(Sp-So)**2.0
            [|deltagreek|]