namespace DerivaRisk.model
open System
open InstrumentData
open DerivaRisk.data
open MathNet.Numerics.Distributions

type MonteCarloGBM(simconfig:SimulationCube,
                   instr:Instrument)=

    inherit MonteCarloModel(simconfig,instr)

    override self.generatepath(md:Marketdata):float[][][]=
        let simdata = self.simconfig
        let (assets,simcube)=MonteCarloCube.generate_cube simdata (int(simdata.seed))
        

        let ntimesteps = int(simdata.number_of_time_steps)
        let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
        let T = float(instr.expiry.Subtract(md.referenceDate).TotalDays)/360.0

        let S = self.optionbase.basket.components
                    |> Array.map(fun id -> id.id,Extractors.extract_spot_for_time_t md t id)
                    |> dict


        let r = Extractors.Integrate_Rate md instr.currency (t,T)
        let q = Extractors.Integrate_Dividends md instr.currency (t,T)
        
        let fvol = assets |> Array.map(fun asset -> Extractors.extract_volatility_surface md asset)

        let optionbase = instr :?> OptionBase
        if optionbase = null then
            failwith "Invalid option type provided to Monte Carlo pricer"

        let strike= optionbase.strike.amount                                  

        //cube is a matrix of correlated npaths x assets x time
        let mc = simcube
                    |> Array.mapi(fun nsim asset_paths ->
                                                       
                            asset_paths
                            
                                
                            )
                
                        )
        
        mc    



    override self.payoff(simpath:float[][]):float=
        failwith "Not Implemented"