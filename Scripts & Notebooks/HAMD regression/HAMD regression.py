import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from scipy.stats import pearsonr
import concurrent.futures
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib
import os

def load_data(filepath):
    df = pd.read_csv(filepath, index_col=0)
    return df.dropna(subset=['HAMD'])

def elasticnet_grid_search(X, y, ratios, alphas, k_splits=5, log_file='elastic_net_cv_progress.txt'):
    kf = KFold(n_splits=k_splits, shuffle=True, random_state=42)
    best_ratio = None
    best_alpha = None
    best_mse = float('inf')
    best_r2 = float('-inf')
    progress = 0
    total_models = len(ratios) * len(alphas) * kf.get_n_splits(X)

    with open(log_file, 'w') as f:
        for ratio in ratios:
            for alpha in alphas:
                for train_idx, val_idx in kf.split(X):
                    X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
                    y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

                    model = ElasticNet(alpha=alpha, l1_ratio=ratio, random_state=42, max_iter=10000)
                    model.fit(X_train, y_train)
                    y_pred = model.predict(X_val)

                    mse = mean_squared_error(y_val, y_pred)
                    r2 = r2_score(y_val, y_pred)

                    if mse < best_mse:
                        best_mse, best_ratio, best_alpha, best_r2 = mse, ratio, alpha, r2

                    progress += 1
                    msg = f"Progress: {progress}/{total_models}, l1_ratio={ratio}, alpha={alpha}, MSE={mse}, R2={r2}"
                    print(msg)
                    f.write(msg + '\n')

    return best_ratio, best_alpha

def bootstrap_feature_selection(X, y, alpha, l1_ratio, n_iterations=100):
    X_array, y_array = X.to_numpy(), y.to_numpy()
    def run_iteration(i):
        idx = np.random.choice(range(X_array.shape[0]), size=int(0.8 * X_array.shape[0]), replace=True)
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=42)
        model.fit(X_array[idx], y_array[idx])
        return model.coef_ != 0

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(run_iteration, range(n_iterations)))

    selection_freq = np.mean(np.array(results).T, axis=1)
    return pd.DataFrame({'Feature': X.columns, 'Frequency': selection_freq})

def find_best_threshold(X, y, feature_freq, alpha, l1_ratio):
    best_threshold = 0
    best_r2 = float('-inf')
    thresholds = np.arange(1.0, 0.0, -0.01)
    for thresh in thresholds:
        mask = feature_freq['Frequency'] >= thresh
        if mask.sum() == 0:
            continue
        X_sel = X.loc[:, mask.values]
        r2s = [r2_score(y_val, ElasticNet(alpha=alpha, l1_ratio=l1_ratio).fit(
                    *train_test_split(X_sel, y, train_size=0.8, random_state=42))
                .predict(X_sel)) for _ in range(10)]
        avg_r2 = np.mean(r2s)
        if avg_r2 > best_r2:
            best_r2, best_threshold = avg_r2, thresh
    return best_threshold

def final_evaluation(X, y, alpha, l1_ratio, threshold, feature_freq):
    mask = feature_freq['Frequency'] >= threshold
    X_selected = X.loc[:, mask.values]

    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    y_true_all, y_pred_all = [], []
    mse_list, r2_list, mae_list, pearson_list = [], [], [], []

    for train_idx, val_idx in kf.split(X_selected):
        X_train, X_val = X_selected.iloc[train_idx], X_selected.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=42)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_val)

        y_true_all.extend(y_val)
        y_pred_all.extend(y_pred)

        mse_list.append(mean_squared_error(y_val, y_pred))
        r2_list.append(r2_score(y_val, y_pred))
        mae_list.append(mean_absolute_error(y_val, y_pred))
        pearson_list.append(pearsonr(y_val, y_pred)[0])

    return np.mean(mse_list), np.mean(r2_list), np.mean(mae_list), np.mean(pearson_list), y_true_all, y_pred_all

def plot_correlation(y_true, y_pred, r2, mae, pearson_corr, save_path):
    df_pred = pd.DataFrame({'True': y_true, 'Pred': y_pred})
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("grass_to_orange", ["#8FBC8F", "#FF8C00"], N=256)
    norm = Normalize(vmin=min(y_true), vmax=max(y_true))
    sm = ScalarMappable(cmap=cmap, norm=norm)
    df_pred['Color'] = df_pred['True'].apply(lambda x: sm.to_rgba(x))

    plt.figure(figsize=(8, 6))
    sns.regplot(x='True', y='Pred', data=df_pred, scatter=False, color='grey')
    sns.boxplot(x='True', y='Pred', data=df_pred, palette=df_pred['Color'])

    plt.xlim(-3, 39)
    plt.xticks(np.arange(-1, 36, 5))
    plt.xlabel('True HAMD score')
    plt.ylabel('Predicted HAMD score')
    plt.title('Correlation')

    plt.text(0.95, 0.25, f'R²: {r2:.3f}', ha='right', va='top', transform=plt.gca().transAxes)
    plt.text(0.95, 0.20, f'MAE: {mae:.3f}', ha='right', va='top', transform=plt.gca().transAxes)
    plt.text(0.95, 0.15, f'Pearson corr: {pearson_corr:.3f}', ha='right', va='top', transform=plt.gca().transAxes)

    plt.savefig(save_path)
    plt.close()

def main():
    df = load_data('./mirna.csv')
    X, y = df.drop(columns=['HAMD']), df['HAMD']

    ratios = np.arange(0.01, 1.01, 0.05)
    alphas = np.concatenate(([10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 0], np.arange(0.1, 1.0, 0.1), np.arange(1, 11, 1)))

    best_ratio, best_alpha = elasticnet_grid_search(X, y, ratios, alphas)
    feature_freq_df = bootstrap_feature_selection(X, y, best_alpha, best_ratio)
    feature_freq_df.to_csv('feature_selection_frequency.csv', index=False)

    best_thresh = find_best_threshold(X, y, feature_freq_df, best_alpha, best_ratio)
    mse, r2, mae, pearson_corr, y_true, y_pred = final_evaluation(X, y, best_alpha, best_ratio, best_thresh, feature_freq_df)

    print(f"\nFinal evaluation with threshold {best_thresh}:")
    print(f"R²: {r2}, MAE: {mae}, Pearson corr: {pearson_corr}")

    plot_correlation(y_true, y_pred, r2, mae, pearson_corr, 'mirna_correlation_distribution.pdf')

if __name__ == '__main__':
    main()
