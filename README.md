# syn_bmca

A repository for hosting the data and code for Synechococcus elongatus using Bayesian Metabolic Control Analysis.

```mermaid
flowchart TD
 subgraph Legend[<h2> Legend </h2>]
 direction LR
        UserData["User Provided Data"]
        CompData["Computed Data"]
        Action(("Action/Process"))
        Choice{"Choice"}
  end

 subgraph ObservedData[<h2> Observed Data and Pre-processing </h2>]
        Data["Observed Data"]
        TranscData["Transcriptomics"]
        ExMetab["External Metabolites"]
        RxnConstr["Reactions of Interest & Objectives"]
        RefSt["Reference State"]
        CalcExFlux["Calculated External Fluxes"]
        EnzAct["Enzyme Activity"]
  end

 subgraph CobraModel[<h2> Cobra Model and Pre-processing </h2>]
        FVA(("FVA"))
        ElastMat["Elasticity Matrices"]
        StoichMat["Stoichiometric Matrix"]
        Optimize(("Optimize"))
        FBA(("FBA"))
        OrigCobra["Cobra Model"]
        FluxBounds["Flux Bounds"]
        RefCobra["Data-Conditional Cobra Model"]
        UncondCompFlux["Unconditional Computed Fluxes"]
        RefStFlux["Reference Strain Fluxes"]
        DelZeroFlux(("Delete 0-Flux Reactions"))
        Eflux(("Eflux2"))
        CondCompFluxes["Conditional Computed Fluxes"]
  end

 subgraph PyMC[<h2> PyMC Model Build </h2>]
        StartPyMC(("Start"))
        ChckObsFlux{"Use Data-Conditioned Fluxes?"}
        GetInputsBasic(("Collect Minimal Inputs"))
        CreatePyMC(("Generate PyMC Model"))
        GetInputsData(("Collect Data-informed Inputs"))
  end

 subgraph Bayes[<h2> Bayesian Analysis </h2>]
        RunInf(("Run Inference"))
        Priors["Prior FCCs"]
        Post["Posterior FCCs"]
        Viz["Visualization and Comparisons"]
  end
    
    %% Data Links
    Data --> TranscData & ExMetab & RxnConstr & RefSt
    ExMetab --> CalcExFlux
    TranscData -- Convert --> EnzAct

    %% Cobra Model Links
    OrigCobra --> FVA & ElastMat & Optimize & FBA
    FVA --> FluxBounds
    FluxBounds --> RefCobra
    Optimize --> UncondCompFlux --> RefStFlux
    FBA --> RefStFlux
    RefSt --> RefStFlux
    RefCobra --> DelZeroFlux & Eflux & StoichMat
    DelZeroFlux --> RefCobra
    Eflux --> CondCompFluxes

    %% PyMC Model Build
    StartPyMC --> ChckObsFlux
    ChckObsFlux -- No --> GetInputsBasic
    StoichMat --> GetInputsBasic & GetInputsData
    ElastMat --> GetInputsBasic & GetInputsData
    RefStFlux --> GetInputsBasic & GetInputsData
    GetInputsBasic --> CreatePyMC
    ChckObsFlux -- Yes --> GetInputsData
    CondCompFluxes --> GetInputsData
    GetInputsData --> CreatePyMC

    %% Inter-subgraph Links
    RxnConstr --> FVA
    CalcExFlux & EnzAct --> Eflux
    CreatePyMC --> RunInf & Priors
    RunInf --> Post
    Priors --> Viz
    Post --> Viz

    %% Color styles for each node
    %% Raw Data - Blue
    style UserData stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style TranscData stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style ExMetab stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style RxnConstr stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style RefSt stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style Data stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style OrigCobra stroke:#2962FF,fill:#2962FF,color:#FFFFFF

    %% Computed Data- Green
    style CompData fill:#00C853,color:#000000
    style CalcExFlux fill:#00C853,color:#000000
    style EnzAct fill:#00C853,color:#000000
    style FVA stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style ElastMat fill:#00C853,color:#000000
    style StoichMat fill:#00C853,color:#000000
    style Priors fill:#00C853,color:#000000
    style Post fill:#00C853,color:#000000
    style Viz fill:#00C853,color:#000000
    style FluxBounds fill:#00C853,color:#000000
    style RefCobra fill:#00C853,color:#000000
    style UncondCompFlux fill:#00C853,color:#000000
    style RefStFlux fill:#00C853,color:#000000
    style CondCompFluxes fill:#00C853,color:#000000

    %% Action/Process - Yellow
    style Action stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style Optimize stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style FBA stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style DelZeroFlux stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style Eflux fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
    style StartPyMC fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
    style GetInputsBasic fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
    style CreatePyMC fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
    style GetInputsData fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
    style RunInf fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000

    %% Choice - Orange
    style Choice stroke:#FF6D00,stroke-width:4px,stroke-dasharray: 0,fill:#FF6D00,color:#000000
    style ChckObsFlux stroke:#FF6D00,stroke-width:4px,stroke-dasharray: 0,fill:#FF6D00,color:#000000
```
