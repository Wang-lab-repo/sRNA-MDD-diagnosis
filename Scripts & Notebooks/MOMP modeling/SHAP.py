import pandas as pd
import matplotlib.pyplot as plt
import shap
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt


output_dir = './'

num_features = X_test.shape[1]
max_evals = 2 * num_features + 1

shap_values_dict = {}

for name, model in best_models.items():
    print(f"Processing SHAP for {name}...")

    if hasattr(model, "predict_proba"):
        explainer = shap.Explainer(model.predict_proba, mirna_feature)
    else:
        explainer = shap.Explainer(model.predict, mirna_feature)

    shap_values = explainer(X_test, max_evals=max_evals)

    if len(shap_values.values.shape) == 3:
        shap_values_selected = shap_values[..., 1]
    else:
        shap_values_selected = shap_values
    shap_values_dict[name] = shap_values_selected
    
    output_path_csv = f'{output_dir}rna_{name}_shap_values.csv'
    pd.DataFrame(shap_values_selected.values, columns=X_test.columns).to_csv(output_path_csv, index=False)
    print(f"SHAP values for {name} saved to {output_path_csv}")


start_color = np.array([0, 51, 102]) 
end_color = np.array([255, 140, 0])

num_colors = 10
gradient_colors = [start_color + (end_color - start_color) * (i / (num_colors - 1)) for i in range(num_colors)]
gradient_colors = np.array(gradient_colors) / 255

cmap = ListedColormap(gradient_colors)

for name in best_models.keys():
    print(f"Generating SHAP plot for {name} from CSV...")

    input_path_csv = f'{output_dir}rna_{name}_shap_values.csv'
    shap_values_selected = pd.read_csv(input_path_csv)

    plt.figure()
    shap.summary_plot(shap_values_selected.values, X_test, max_display=5, show=False, cmap = cmap)

    output_path_pdf = f'{output_dir}rna_{name}_shap.pdf'
    plt.savefig(output_path_pdf, format='pdf')
    plt.close()
    print(f"SHAP plot for {name} saved to {output_path_pdf}")