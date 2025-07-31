from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize
import numpy as np

class ModelSelectionProblem(ElementwiseProblem):
    def __init__(self, data_dict, models, objectives):
        self.data_dict = data_dict
        self.models = models
        self.objectives = objectives
        super().__init__(n_var=len(models), n_obj=len(objectives), xl=0.0, xu=1.0, type_var=float)

    def _evaluate(self, x, out, *args, **kwargs):
        x_binary = np.clip(np.round(x).astype(int), 0, 1)
        selected_models = [self.models[i] for i, val in enumerate(x_binary) if val == 1]

        if not selected_models:
            out["F"] = [1.0] * len(self.objectives)
            return

        scores = []
        for obj in self.objectives:
            obj_scores = []
            for dataset, data in self.data_dict.items():
                selected_data = data[data["Model"].isin(selected_models)]
                if selected_data.empty:
                    obj_scores.append(1.0)
                else:
                    mean_val = selected_data[obj].mean()
                    if obj in ["Log-Loss", "InferenceTime"]:
                        obj_scores.append(mean_val)
                    else:
                        obj_scores.append(-mean_val)
            scores.append(np.mean(obj_scores))
        out["F"] = scores

def optimize_model_selection(performance_df, objectives):
    from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting

    data_dict = {
        ds: performance_df[performance_df["Dataset"] == ds]
        for ds in performance_df["Dataset"].unique()
    }

    models = performance_df["Model"].unique().tolist()
    problem = ModelSelectionProblem(data_dict, models, objectives)

    ref_dirs = get_reference_directions("das-dennis", n_dim=len(objectives), n_partitions=6)
    algorithm = NSGA3(pop_size=100, ref_dirs=ref_dirs)

    res = minimize(problem, algorithm, ('n_gen', 100), seed=42, verbose=True)

    nds = NonDominatedSorting().do(res.F, only_non_dominated_front=True)
    pareto_X = res.X[nds]
    pareto_F = res.F[nds]

    pareto_records = []
    for i, (x, f) in enumerate(zip(pareto_X, pareto_F)):
        x_binary = np.round(x).astype(int)
        selected = [models[i] for i, val in enumerate(x_binary) if val == 1]
        record = {
            "Selected_Models": selected,
        }
        for obj_name, obj_val in zip(objectives, f):
            if obj_name in ["Log-Loss", "InferenceTime"]:
                record[obj_name] = obj_val
            else:
                record[obj_name] = -obj_val
        pareto_records.append(record)

    df_pareto = pd.DataFrame(pareto_records)
    df_pareto.to_csv('./result/pareto_model_pool.csv', index=False)
    print(f"\nSaved {len(df_pareto)} Pareto-optimal solutions to ./result/pareto_model_pool.csv")

    return df_pareto