(* ::Package:: *)

(* ::Text:: *)
(*Chemical Reaction Network (CRN) Simulator package is developed by David Soloveichik. Copyright 2009-2015. *)
(*http://users.ece.utexas.edu/~soloveichik/crnsimulator.html*)


Needs["CRNSimulator`"];
BeginPackage["CRNSimulatorSSA5`", {"CRNSimulator`"}];

(* Warning:  species names like A[3][17] are not supported *)

(* with modifications by Erik Winfree, 2016 *)
(* UPDATE v2: SSAToListPlot has been rewritten to be time-based rather than event-based. *)
(* UPDATE v2: computePropensities has been rewritten to compile better, and now uses rxnsysToSimulationStructures. *)
(* UPDATE v2: SimulateRxnsysSSA now has options to either return a full trajectory or just the first/last states. *)
(* UPDATE v3: rxnsysToSimulationStructures now works out which other reactions' propensities are affected by each event. *)
(* UPDATE v3: computePropensities now just does one specified reaction, and is used to only update affected reactions. *)
(* UPDATE v4: uses a fixed-array binary tree for propensities to get O(log N) move selections. *)
(* 
    For 50,000 reactions, 5000 species: 
    v3 is 20% slower than v2 for pre-processing, but 20x faster running
    v4 is the same as v3 for pre-processing, and 2x faster than v3.
*)
(* UPDATE v5: allows CompileRxnsysSSA to do the precomputing once, and then run it many times with different initial conditions. *)
(*     ??ALSO?? should I add catrxn, statevector->concs, and other useful functions? some of these should go to CRNSimulator.m??  *)

(* See "Efficient Formulations for Exact Stochastic Simulation of Chemical Systems" Mauch & Stalzer 2011 *)
(* See "Adaptive tree-based search for Stochastic simulation algorithm" Thanh & Zunino 2014 *)


(* ::Section:: *)
(*Public interface specification*)


SimulateRxnsysSSA::usage="SimulateRxnsysSSA[rxnsys,endtime] runs Gillespie's SSA simulation of the reaction \
system rxnsys for time 0 to endtime or until no further reactions are possible. Set endtime to \[Infinity] to continue \
until no further reactions are possible. Alternatively, {begintime, endtime} may be used to specify a non-zero initial time, \
or {begintime, endtime, timeinterval} to collect a trajectory containing the system state sampled at regular intervals (plus the final state). \
In rxnsys, reactions are specified by rxn statements (arbitrary \
order reactions are allowed: eg rxn[a+a+b,c+c,1], including zero-order reactions), and initial molecular counts \
are specified by conc statements. \
If no initial count is set for a species, its initial count is set to 0. \
The shorthand revrxn[x,a+b,3.7,5.2] expands to two irreversible reactions. \
Compatible conc statements and rxn statements are summed. \
SimulateRxnsysSSA returns {times, data} where times is a list of reaction firing times, \
and data is the corresponding list of states. \
The propensity of reaction a\[Rule]\[Ellipsis] is k*a, the propensity of reaction a+b\[Rule]\[Ellipsis] is k*a*b, \
and the propensity of reaction a+a\[Rule]\[Ellipsis] is k*a*(a-1), etc.\
A zero-order reaction can be written e.g. rxn[0,x+y,1]. \
Options: NumStepsLimit \[Rule] _Integer can be used to reduce the memory overhead per batch of simulation steps. \
Output \[Rule] \"State\" stores and returns only the first and last states, thus also reducing memory overhead. \
Output \[Rule] \"Trajectory\" is the default. ";

SSAToListPlot::usage="SSAToListPlot[{times,data}] takes the output of SimulateRxnsysSSA and converts it \
to a form that can be directly plotted using ListPlot, ListLinePlot, etc. \
SSAToListPlot[{times,data},maxpoints] limits the number of points to plot to maxpoints by sampling \
at evenly-spaced time intervals (useful for long runs).";

CompiledRxnsysSSA::usage="CompiledRxnsysSSA[rxnsys] precomputes data structures needed by the simulator, \
so that the same system can be efficiently run many times with different initial counts. \
After precomputation, the format is CompiledRxnsysSSA[rxnsys,datastructures].  This format can be given to \
SimulateRxnsysSSA in the place of a standard-format rxnsys.";

ConcentrationsFromSSAState::usage="ConcentrationsFromSSAState[rsys,state] converts a counts vector \
to a list of conc[species,count] entries, ignoring zero-count species.";

NewConcentrations::usage="NewConcentrations[rxnsys,newconcs] efficiently changes the concentrations \
for a compiled or non-compiled reaction system.";



(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* Set Compiler options *)
On["CompilerWarnings"];
System`SetSystemOptions["CompileOptions"->"CompileReportExternal"->True];
On[Compile::noinfo];
SetOptions[Compile,CompilationOptions->{"InlineExternalDefinitions" -> True,"InlineCompiledFunctions"->False},
CompilationTarget->"C",Parallelization->True,RuntimeOptions->"Speed"];


(* Produces reactions in srxn[[reaction index],{x1,x2},{x3},k] format ("structured reaction system"). Removes any other statements. *)
rxnsysToSrxnsys[rsys_] := 
	Module[{i = 1}, 
		Cases[
			rsys, 
			rxn[rs_, ps_, k_] :> 
				srxn[i++, 
					{Unevaluated[rs] /. {Plus -> Seq, c_Integer*s_ :> Seq @@ Table[s, {c}]}}, {Unevaluated[ps] /. {Plus -> Seq, c_Integer*s_ :> Seq @@ Table[s, {c}]}}, 
					k]]]


rxnsysToSimulationStructures[rsys_] := Module[{spcs, hash, Nrxns,Nspcs,maxRA,maxPA,maxCA,maxNA, srxnsys,
      symbolicReactants, symbolicProducts, symbolicChange, species, indices, stoichs,
      affectedby, reactionsAffectedTemp,
      rateConstants,
	  productionArities, reactionArities, reactionSpecies, reactionStoichs,
	  stateChangeVectors, stateChangeArities, stateChangeSpecies, stateChangeStoichs,
	  numAffected, reactionsAffected},
    spcs = SpeciesInRxnsys[rsys]; Nspcs = Length[spcs];
    Do[hash[spcs[[i]]]=i,{i,Nspcs}];
	srxnsys=rxnsysToSrxnsys[rsys];
    rateConstants=Cases[srxnsys,srxn[_,_,_,k_]:>k];
    symbolicReactants = Cases[srxnsys, srxn[_, r_, p_, _] :> Plus @@ r];
    symbolicProducts = Cases[srxnsys, srxn[_, r_, p_, _] :> Plus @@ p];
	symbolicChange = Cases[srxnsys, srxn[_, r_, p_, _] :> (Plus @@ p)-(Plus @@ r)];
    species[sr_]:= Replace[MonomialList[sr], {Times[a_,b_]:>b,_Integer->Seq[]},1];
    indices[sr_]:=(hash /@ species[sr]);
    stoichs[sr_]:=(Coefficient[sr,#]& /@ species[sr]);
    reactionArities = Length[species[#]]& /@ symbolicReactants;
	productionArities = Length[species[#]]& /@ symbolicProducts;
    Nrxns = Length[reactionArities]; maxRA = Max[reactionArities]; maxPA = Max[productionArities];
    reactionSpecies = Table[0, {n, Nrxns}, {i, maxRA}];
    Do[reactionSpecies[[n, Range[reactionArities[[n]]]]] = indices[symbolicReactants[[n]]], {n, Nrxns}];
    reactionStoichs = Table[0, {n, Nrxns}, {i, maxRA}];
    Do[reactionStoichs[[n, Range[reactionArities[[n]]]]] = stoichs[symbolicReactants[[n]]], {n, Nrxns}];
(*
Print["rateConstants = ",rateConstants];    
Print["reactionArities = ",reactionArities];
Print["reactionSpecies = ",reactionSpecies];
Print["reactionStoichs = ",reactionStoichs];
*)
    stateChangeVectors = Table[0,{n,Nrxns},{i, Nspcs}];
    Do[stateChangeVectors[[n,indices[symbolicReactants[[n]]]]] -= stoichs[symbolicReactants[[n]]],{n,Nrxns}];
    Do[stateChangeVectors[[n,indices[symbolicProducts[[n]]]]] += stoichs[symbolicProducts[[n]]],{n,Nrxns}];

	stateChangeArities = Length[species[#]]& /@ symbolicChange; maxCA = Max[stateChangeArities];
    stateChangeSpecies = Table[0, {n, Nrxns}, {i, maxCA}];
    Do[stateChangeSpecies[[n, Range[stateChangeArities[[n]]]]] = indices[symbolicChange[[n]]], {n, Nrxns}];
    stateChangeStoichs = Table[0, {n, Nrxns}, {i, maxCA}];
    Do[stateChangeStoichs[[n, Range[stateChangeArities[[n]]]]] = stoichs[symbolicChange[[n]]], {n, Nrxns}];
(*
Print["stateChangeVectors = ", stateChangeVectors];
Print["stateChangeArities = ", stateChangeArities];
Print["stateChangeSpecies = ", stateChangeSpecies];
Print["stateChangeStoichs = ", stateChangeStoichs];
*)
	(* which reactions' propensities will change if species s changes? *)
	Do[affectedby[s]={},{s,Nspcs}]; 
    Do[AppendTo[affectedby[#],n]& /@ reactionSpecies[[n,Range[reactionArities[[n]]] ]], {n,Nrxns}];

	(* which reactions need their propensities updated after reaction n has occurred? *)
	reactionsAffectedTemp = Table[
		Union @@ (affectedby /@ stateChangeSpecies[[n,Range[stateChangeArities[[n]]]]]),
	    {n,Nrxns}];
	(* compiled Mathematica functions only accept rectangular arrays *)
	numAffected = Length /@ reactionsAffectedTemp;  maxNA = Max[numAffected];
(*
Print["affectedby computed: ", spcs, " -> ", affectedby /@ Range[Nspcs]];
Print["reactionsAffectedTemp = ", reactionsAffectedTemp];
Print["numAffected = ", numAffected];
Print["maxNA = ", maxNA];
*)
	reactionsAffected = Table[0, {n, Nrxns}, {i, maxNA}];
    Do[reactionsAffected[[n, Range[numAffected[[n]]]]] = reactionsAffectedTemp[[n]], {n, Nrxns}];

    {stateChangeVectors,1.0*rateConstants,
     reactionArities,reactionSpecies,reactionStoichs,
     numAffected, reactionsAffected}
]


computeWays = Compile[{{count,_Integer},{stoich,_Integer}}, If[count<stoich,0,Product[c, {c,count-stoich+1,count} ]]];

computePropensities = Compile[{{w,_Integer,1},{rxn,_Integer},{rates,_Real,1},{arity,_Integer,1},{species,_Integer,2},{stoichs,_Integer,2}},
           rates[[rxn]]*Product[ computeWays[w[[species[[rxn,i]]]],stoichs[[rxn,i]]],{i,arity[[rxn]]}]
         ];

SSAengine = Compile[{{stateChangeVectors,_Integer,2}, {Naffected,_Integer,1}, {affected,_Integer,2},
       {rates,_Real,1}, {arity,_Integer,1}, {species,_Integer,2}, {stoichs,_Integer,2},
       {x0, _Integer, 1}, {t0, _Real}, {tlimit, _Real}, 
       {expRvPool, _Real, 1}, {unifRvPool, _Real, 1},
       {outState,True|False},{outTraj,True|False}},
     
     Module[{x = x0, t = t0, propensities, sumpropensities, treeDepth,
       rvPointer = 1, rxnNum = Length[stateChangeVectors], 
       spcsNum = Length[First[stateChangeVectors]], nextRxn, 
       maxNumSteps = Length[expRvPool], deltat, output, outsize = 2, step = 1, s, 
       mark, d, j},
 
      (* the first dimension is time *)
      If[outTraj, outsize=maxNumSteps];  
      If[outState, outsize=2];  (* same format as output, but just 1 entry (plus extra) *)
      output = Table[0.0, {outsize}, {spcsNum + 1}];
      
	  (* Each node in the tree is the sum of the propensities of its children, i.e.,
          propensities[[treeDepth-d,i]] = 
             Sum[ propensities[[treeDepth-d+1,j]], {j,(i-1)*2^d + 1, (i)*2^d}] *)
	  treeDepth = Max[1,Ceiling[Log2[1.0*rxnNum]]]; (* compiler does typecast _Int \[Rule] _Real here! *)
	  propensities = Table[0.0, {d, treeDepth},{rxn,2^treeDepth}];
(*
Print["rxnNum = ", rxnNum, ", treeDepth = ", treeDepth,", propensities sizes = ", Length /@ propensities];
Print["Rates = ", rates];
propensities[[treeDepth,5]] = computePropensities[x,5,rates,arity,species,stoichs];
Print["End Test"];
*)
	  Do[propensities[[treeDepth,rxn]] = computePropensities[x,rxn,rates,arity,species,stoichs],
              {rxn,rxnNum}];
(*
Print["Reaction propensities initialized for ", rxnNum, " reactions, x = ", x];
*)
	  Do[ Do[ propensities[[d,i]] = propensities[[d+1,2*i-1]]+propensities[[d+1,2*i]], 
              {i,2^d}], {d,treeDepth-1,1,-1}];

      While[
       (* use tlimit or maxNumSteps only if nonzero *)
       If[tlimit == 0, True, t <= tlimit] && 
        If[maxNumSteps == 0, True, step <= maxNumSteps], 
       
       sumpropensities = propensities[[1,1]]+propensities[[1,2]];
       If[sumpropensities == 0.0, Break[]];  (* no reaction is possible *)
       
       (* compute deltat and nextRxn*)
       
       deltat = (1/sumpropensities)*expRvPool[[rvPointer]];
       mark = sumpropensities*unifRvPool[[rvPointer]];
       rvPointer++;

	   (* walk down the binary tree to choose the move *)
       For[j=1; d=1, d<=treeDepth, d++,
           If[mark < propensities[[d,j]], 
              j=2*j-1,
              mark -= propensities[[d,j]]; j=2*j+1
           ]];
       nextRxn = Min[Max[1,Quotient[j+1,2]],rxnNum];

       (* implement reaction *)
       
       x += stateChangeVectors[[nextRxn]];
       t += deltat;

	   (* update propensities for the reactions that are affected *)
	   Do[ propensities[[treeDepth,affected[[nextRxn,n]]]] = 
			    computePropensities[x,affected[[nextRxn,n]],rates,arity,species,stoichs],
		   {n, Naffected[[nextRxn]]}];
       (* maintain consistency of the move tree *)
	   Do[ For[j=affected[[nextRxn,n]]; d=treeDepth, d>1, d--,
			 j = Quotient[j+1,2];
             propensities[[d-1,j]] = propensities[[d,2*j-1]]+propensities[[d,2*j]] ],
		   {n, Naffected[[nextRxn]]}];
(*
Print["t = ",t, " : rxn = ", nextRxn, " : x = ", x, " : prop = ", propensities[[treeDepth]]];
Do[Print[propensities[[d,Range[2^d]]]],{d,treeDepth}];
*)
       
	   If[outTraj,
         output[[step, 1]] = t;
         For[s = 1, s <= spcsNum, s++, output[[step, 1 + s]] = x[[s]] ]
       ];
       
       step++;
      ];
      
      If[outState, 
		output[[1,1]]=t;
        For[s = 1, s <= spcsNum, s++, output[[1, 1 + s]] = x[[s]] ];
        output[[2,1]]= step-1;  (* number of events simulated *)
        step=3; (* so we output 2 data points, one being the above *)
      ];
      output[[;; step - 1]]
      ]
    ];

CompiledRxnsysSSA[rsys_List] := CompiledRxnsysSSA[rsys,rxnsysToSimulationStructures[rsys]]

ConcentrationsFromSSAState[rsys_List,data_List] := Cases[MapThread[conc[#1, #2] &, {SpeciesInRxnsys[rsys], Last[data]}], conc[_,_?Positive]]
ConcentrationsFromSSAState[CompiledRxnsysSSA[rsys_List,datastructures_List],data_List] := ConcentrationsFromSSAState[rsys,data]

NewConcentrations[rsys_List,concs_List] := Join[Cases[rsys,_rxn],concs]
NewConcentrations[CompiledRxnsysSSA[rsys_List,datastructures_List],concs_List] := CompiledRxnsysSSA[NewConcentrations[rsys,concs],datastructures]

SpeciesInRxnsys[CompiledRxnsysSSA[rsys_List,datastructures_List]] := SpeciesInRxnsys[rsys]

Options[SimulateRxnsysSSA] = {NumStepsLimit->0,Output->"Trajectory"};

SimulateRxnsysSSA[rsysOrCompiled_, tlimitarg_, OptionsPattern[]] := 
 Module[{rsys, concs, spcs, x0arg, begintime=0.0, endtime=Infinity, timeinterval=0.0,
    stateChangeVectors, rateConstants, reactionArities, reactionSpecies, reactionStoichs, numAffected, reactionsAffected, 
    outputTypes, outputTrajectory=False, outputState=True, 
    w, i, stepsSimulated, totalSteps=0, t0, t1,
    maxNumStepsarg, expRvPool, unifRvPool, 
    out, times, data, timesrun, datarun},
  
  If[Head[rsysOrCompiled]===CompiledRxnsysSSA,
     rsys=First[rsysOrCompiled];
     {stateChangeVectors,rateConstants,reactionArities,reactionSpecies,reactionStoichs,numAffected,reactionsAffected} =
          Last[rsysOrCompiled],
     (* else *)
	 rsys=rsysOrCompiled;
     t0=TimeUsed[];
     {stateChangeVectors,rateConstants,reactionArities,reactionSpecies,reactionStoichs,numAffected,reactionsAffected}=
          rxnsysToSimulationStructures[rsys];
     t1=TimeUsed[];
     PrintTemporary["Took ",t1-t0," seconds to preprocess the ",Length[rateConstants],"-reaction CRN."];
     PrintTemporary["Preprocessing done, starting simulation..."];
  ];
  PrintTemporary["Maximum number of reactions affected by another reaction's execution: ",Max[numAffected]];

  spcs = SpeciesInRxnsys[rsys];

  outputTypes = OptionValue[Output];
  outputTrajectory=(outputTypes=="Trajectory" || MemberQ[outputTypes,"Trajectory"]);
  outputState=(outputTypes=="State" || MemberQ[outputTypes,"State"]);
  maxNumStepsarg = With[{m=OptionValue[NumStepsLimit]}, 
     If[m==0, If[outputTrajectory,Ceiling[10^6/Length[spcs]],10^6], m]];

  (* Vector of initial counts from parsing conc statements.  Multiple conc for same species are summed. *)
  x0arg = Plus @@ Cases[rsys, conc[#, c_] :> c] & /@ spcs;
  
  (* initial state as first point *)
  Switch[ Length[tlimitarg],
    0, endtime=tlimitarg, (* integer or real value *)
	1, endtime=tlimitarg, (* Infinity *)
    2, begintime=First[tlimitarg]; endtime=Last[tlimitarg],
    3, begintime=First[tlimitarg]; endtime=tlimitarg[[2]]; timeinterval=Last[tlimitarg]; Print["WARNING: time intervals not yet implemented."],
    _, Print["WARNING: time limit argument unrecognized, simulating to Infinity."]
  ];
  times = {begintime};
  data = {x0arg};
  
  While[True,
   expRvPool = RandomVariate[ExponentialDistribution[1], maxNumStepsarg];
   unifRvPool = RandomVariate[UniformDistribution[], maxNumStepsarg];
   out = SSAengine[stateChangeVectors, numAffected, reactionsAffected, 
     rateConstants, reactionArities, reactionSpecies, reactionStoichs,
     Last[data], Last[times], endtime /. \[Infinity] -> 0, 
     expRvPool, unifRvPool,
     outputState, outputTrajectory];
   If[outputTrajectory,
     timesrun = out[[All, 1]];
     datarun = Round[out[[All, 2 ;;]]];  
     times = Join[times, timesrun];
     data = Join[data, datarun];
	 stepsSimulated = Length[timesrun]
   ];
   If[outputState,  (* keep it a trajectory of length 2 *)
     times = {0.0, out[[1,1]]};
     data = Round[{x0arg, out[[1,2;;]]}];
     stepsSimulated = Round[out[[2,1]]]
   ];
   If[stepsSimulated<maxNumStepsarg, Break[]];
   totalSteps+=stepsSimulated;
   PrintTemporary[totalSteps," events simulated, t = ", Last[times]," ..."];
  ];

  {times, data} 
  
]


(* EW edits this to change from evenly-spaced events to evenly-spaced time points *)
SSAToListPlot[{times_,data_},maxpoints_:\[Infinity]]:=
	Module[{deltaplot,timesplot,dataplot,epochs,epochal,indices},

	(* Decrease the number of points to plot to at most maxpoints*)
	If[Length[times] <= maxpoints,
		timesplot=times; dataplot=data,
		deltaplot = (Last[times]-First[times])/(maxpoints-1);
		epochs = Floor[(times-First[times])/deltaplot];
		epochs[[-1]]=maxpoints-1; (* corrects potential round-off error *)
        epochal = Differences[epochs]; 
		indices = Join[{1},Flatten[MapIndexed[Table[First[#2],#1]&,epochal]]];
		timesplot = Table[First[times]+n*deltaplot,{n,0,maxpoints-1}];
		dataplot = data[[indices]];
        ];

	Transpose[{timesplot,#}]&/@Transpose[dataplot]]
    

(* Reset Compile options so Plot, etc works ok *)
SetOptions[Compile,CompilationOptions->All];


End[];
EndPackage[];
