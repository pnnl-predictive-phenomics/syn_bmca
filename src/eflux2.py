import numpy as np
from optlang.symbolics import add

def EFlux2(model, Transcriptomics):
    eflux2_model = model.copy()
    # Parse GPR into a dict containing isozymes (separated by 'or')
    # Each isozyme has a set of subunits (separated by 'and')
    gpr_dict = dict()
    for r in eflux2_model.reactions:
        if r.gene_reaction_rule:
            temp = set()
            for x in [x.strip('() ') for x in r.gene_reaction_rule.split(' or ')]:
                temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
            gpr_dict[r.id] = temp
    # Set the bounds using the transcriptomics data
    for r in eflux2_model.reactions:
        if r.gene_reaction_rule:
            t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
                                else np.array([np.Inf]) for g in p])
                        for p in gpr_dict[r.id]])
            if r.lower_bound < 0.0:
                r.lower_bound = -t
            else:
                pass
            if r.upper_bound > 0.0:
                r.upper_bound = t
            else:
                pass
        else:
            if r.lower_bound <= -1000.0:
                r.lower_bound = -np.Inf
            if r.upper_bound >= 1000.0:
                r.upper_bound = np.Inf
    # solve FBA to calculate the maximum biomass
    eflux2_model.tolerance = 1e-9
    fba_sol = eflux2_model.optimize()
    print('FBA status', fba_sol.status)
    print('FBA solution', fba_sol.objective_value)
    # Constrain the biomass to the optimal value
    for r in eflux2_model.reactions:
        if r.objective_coefficient:
            r.lower_bound = fba_sol.objective_value
    # minimize the sum of squared flux values
    eflux2_model.objective = eflux2_model.problem.Objective(add([r.flux_expression**2 for r in eflux2_model.reactions]), direction='min')
    eflux2_sol = eflux2_model.optimize()
    print('EFlux2 status', eflux2_sol.status)
    print('EFlux2 solution', eflux2_sol.objective_value)
    # return eflux2 solution
    return(eflux2_sol)