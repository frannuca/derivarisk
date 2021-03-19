namespace DerivaRisk.model
open System
open InstrumentData
open DerivaRisk.data
open MathNet.Numerics.Distributions

type Kparams={K0:float;K1:float;K2:float;K3:float;K4:float;A:float}

type HestonParams={
    S0:float;
    V0:float;
    kappa:float;
    lambda:float;        
    rate:float;
    dividends:float;
    rho:float;
    sigma:float;
    theta:float;
    dt:float;
    gamma_1:float;
    gamma_2:float;
    threshold_phi:float
} with member self.E = exp(-self.kappa*self.dt)


type MonteCarloHeston(simconfig:SimulationCube,instr:Instrument)=
    inherit MonteCarloModel(simconfig,instr)
    
    let computeK(v:HestonParams):Kparams=
           let K0 = -v.kappa*v.rho*v.theta/v.sigma*v.dt
           let K1 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_1*v.dt-v.rho/v.sigma
           let K2 = (v.kappa*v.rho/v.sigma-0.5)*v.gamma_2*v.dt+v.rho/v.sigma
           let K3 = (1.0-v.rho**2.0)*v.gamma_1*v.dt
           let K4 = (1.0-v.rho**2.0)*v.gamma_2*v.dt
           let A = K2+0.5*K4
           {K0=K0;K1=K1;K2=K2;K3=K3;K4=K4;A=A}

    
    override self.generatepath(md:Marketdata)=
        let simdata = self.simconfig
        let (assets,simcube)=MonteCarloCube.generate_cube simdata (int(simdata.seed))
        let (_,simcube_vol) = MonteCarloCube.generate_cube simdata (int(simdata.seed2))

        let ntimesteps = int(simdata.number_of_time_steps)
        let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
        let T = float(instr.expiry.Subtract(md.referenceDate).TotalDays)/360.0

        let S = self.optionbase.basket.components |> Array.map(fun id -> id.id,Extractors.extract_spot_for_time_t md t id) |> dict

        let r = Extractors.Integrate_Rate md instr.currency (t,T)
        let q = Extractors.Integrate_Dividends md instr.currency (t,T)
        
        let hestonparams = self.optionbase.basket.components
                            |> Array.map(fun id -> id.id,Extractors.extract_heston md id)
                            |> Array.map(fun (idx,p) ->idx, {S0=S.[idx];
                                                             V0=p.Vo;
                                                             kappa=p.kappa;
                                                             lambda=p.lambda;
                                                             rate=r;
                                                             dividends=q;
                                                             rho=p.rho;
                                                             sigma=p.sigma;
                                                             theta=p.theta;
                                                             dt=float simconfig.dt;
                                                             gamma_1=p.gamma1;
                                                             gamma_2=p.gamma2;
                                                             threshold_phi=1.5})
                            |>Map.ofArray

        let Ks = hestonparams |> Map.map(fun k v -> computeK v)

        let optionbase = instr :?> OptionBase
        if optionbase = null then
            failwith "Invalid option type provided to Monte Carlo pricer"

        let strike= optionbase.strike.amount

        let rec compute_asset_path(x:HestonParams,K:Kparams,V:float[],S:float[],dZ_asset:float[],dZ_vol:float[])(n:int):float[]=
            if(n>=S.Length) then
                S
            else                
                let dZs= dZ_asset.[n]
                let dZv = dZ_vol.[n]
                let m=x.theta + (V.[n-1]-x.theta)*x.E
                let m2 = m**2.0
                let s2 = V.[n-1]*(x.sigma**2.0)*x.E/x.kappa*(1.0-x.E)+ x.theta*(x.sigma**2.0)/(2.0*x.kappa)*(1.0-x.E)**2.0
                let phi = s2/m2
                let K0 =
                        if phi <= x.threshold_phi then
                            let b= Math.Sqrt(2.0/phi-1.0+Math.Sqrt(2.0/phi*(2.0/phi-1.0)))
                            let a= m/(1.0+b**2.0)
                            
                            V.[n] <- a*(b+dZv)**2.0
                            if K.A<1.0/(2.0*a) then
                                let mx=Math.Exp(K.A*b**2.0*a/(1.0-2.0*K.A*a))/Math.Sqrt(1.0-2.0*K.A*a)
                                -Math.Log(mx)-(K.K1+0.5*K.K3)*V.[n-1]
                            else
                                K.K0
                        else
                            let Uv = Normal.CDF(0.0,1.0,dZv)
                            let p = (phi - 1.0)/(phi+1.0)
                            let beta = (1.0-p)/m
                            V.[n] <-
                                    if Uv <= p then
                                        0.0
                                    else
                                        1.0/beta*Math.Log((1.0-p)/(1.0-Uv))

                            if K.A<beta then
                                let mx=p+beta*(1.0-p)/(beta-K.A)
                                -Math.Log(mx)-(K.K1+0.5*K.K3)*V.[n-1]
                            else
                                K.K0
                
                S.[n]<-S.[n-1]*Math.Exp((x.rate-x.dividends)*x.dt+K0+K.K1*V.[n-1]+K.K2*V.[n]+Math.Sqrt(K.K3*V.[n-1]+K.K4*V.[n])*dZs)
                compute_asset_path(x,K,V,S,dZ_asset,dZ_vol)(n+1)                                       

        //cube is a matrix of correlated npaths x assets x time
        let mc = simcube
                    |> Array.mapi(fun nsim asset_paths ->
                                                       
                            asset_paths
                            |> Array.mapi(fun nasset s_path ->
                                let asset = assets.[nasset]
                                let heston_par = hestonparams.[asset]
                                let Kparam = Ks.[asset]
                                let V = Array.zeroCreate(ntimesteps)
                                let S = Array.zeroCreate(ntimesteps)
                                V.[0]<-heston_par.V0
                                S.[0]<-heston_par.S0                                
                                compute_asset_path(heston_par,Kparam,V,S,s_path,simcube_vol.[nsim].[nasset])(1)
                                
                            )
                
                        )
        
        mc    
                              

    override self.payoff(cube:float[][]):float=
       
        let assets = simconfig.assets_correlation.assets
        let process_components(aggregator:float[]->float)(arr:float[][])=
            [|0 .. arr.[0].Length-1|]
            |> Array.map(fun n -> [|0 .. arr.Length-1|] |> Array.map(fun iasset -> arr.[iasset].[n]) |> aggregator)

        let basketpath=
            match self.optionbase.basket.basketAggregation with
            |BasketAggregation.WORSTOF ->  cube |> process_components Array.min
            |BasketAggregation.BESTOF ->   cube |> process_components Array.max
            |BasketAggregation.AVERAGE ->  cube |> process_components Array.average
            |_ -> failwith "Invalid Basket Aggregation Type"


        let ispayoffAlive= if self.optionbase.barrier <> null then
                            self.optionbase.barrier                             
                                |> Array.map(fun b -> match b.barrierType with
                                                        |BarrierType.UpIn -> basketpath |> Array.exists(fun s -> s>=b.level.amount)
                                                        |BarrierType.DownIn -> basketpath |> Array.exists(fun s -> s<=b.level.amount)
                                                        |BarrierType.UpOut -> basketpath |> Array.forall(fun s -> s<b.level.amount)
                                                        |BarrierType.DownOut -> basketpath |> Array.forall(fun s -> s>b.level.amount)
                                                        |_ -> failwith "Invalid Barrier Type"
                                ) |> Array.fold(fun acc s -> acc && s) true
                            else
                                true

        let ST = basketpath.[basketpath.Length-1]

        match (ispayoffAlive,self.optionbase.optionType) with
        |(false,_)->0.0
        |(_,OptionType.CALL) ->  max (ST - self.optionbase.strike.amount) (0.0)
        |(_,OptionType.PUT) ->  max (self.optionbase.strike.amount-ST) (0.0)
        |_ -> failwith "Option Type not implemented"
        