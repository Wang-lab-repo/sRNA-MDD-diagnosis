import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from sklearn.linear_model import ElasticNet
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from scipy.stats import pearsonr

import concurrent.futures
import os

def load_data(filepath):
    df = pd.read_csv(filepath, index_col=0)
    return df.dropna(subset=['HAMD'])

def elasticnet_grid_search(X, y, ratios, alphas, k_splits=5, log_file='elastic_net_cv_progress.txt'):
    kf = KFold(n_splits=k_splits, shuffle=True, random_state=42)
    best_ratio, best_alpha = None, None
    best_mse = float('inf')

    with open(log_file, 'w') as f:
        for ratio in ratios:
            for alpha in alphas:
                mse_fold = []
                for train_idx, val_idx in kf.split(X):
                    X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
                    y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
                    model = ElasticNet(alpha=alpha, l1_ratio=ratio, random_state=42, max_iter=10000)
                    model.fit(X_train, y_train)
                    y_pred = model.predict(X_val)
                    mse_fold.append(mean_squared_error(y_val, y_pred))
                avg_mse = np.mean(mse_fold)
                msg = f"l1_ratio={ratio}, alpha={alpha}, avg_MSE={avg_mse}"
                print(msg)
                f.write(msg + '\n')
                if avg_mse < best_mse:
                    best_mse = avg_mse
                    best_ratio, best_alpha = ratio, alpha
    return best_ratio, best_alpha

def bootstrap_feature_selection(X, y, alpha, l1_ratio, n_iterations=100):
    def run_iteration(seed):
        idx = np.random.choice(range(X.shape[0]), size=int(0.8 * X.shape[0]), replace=True)
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=seed)
        model.fit(X.iloc[idx], y.iloc[idx])
        return model.coef_ != 0

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(run_iteration, range(n_iterations)))
    
    selection_freq = np.mean(np.array(results), axis=0)
    return pd.DataFrame({'Feature': X.columns, 'Frequency': selection_freq})

def find_best_threshold(X, y, feature_freq_df, alpha, l1_ratio, n_trials=101):
    thresholds = np.arange(1.0, 0.0, -0.01)
    best_threshold = None
    best_r2 = float('-inf')

    for threshold in thresholds:
        selected_features = feature_freq_df[feature_freq_df['Frequency'] >= threshold]['Feature']
        if selected_features.empty:
            print(f"Warning: No features selected at threshold {threshold}. Skipping...")
            continue

        X_selected = X[selected_features]
        r2s, mses = [], []

        for i in range(n_trials):
            X_train, X_val, y_train, y_val = train_test_split(X_selected, y, train_size=0.8, random_state=i)
            model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_val)

            r2s.append(r2_score(y_val, y_pred))
            mses.append(mean_squared_error(y_val, y_pred))

        avg_r2, avg_mse = np.mean(r2s), np.mean(mses)
        print(f"Threshold: {threshold:.2f}, Avg R²: {avg_r2:.4f}, Avg MSE: {avg_mse:.4f}")
        if avg_r2 > best_r2:
            best_r2 = avg_r2
            best_threshold = threshold
    return best_threshold

def final_evaluation(X, y, alpha, l1_ratio, threshold, feature_freq_df):
    selected_features = feature_freq_df[feature_freq_df['Frequency'] >= threshold]['Feature']
    X_selected = X[selected_features]

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
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom", ["#8FBC8F", "#FF8C00"])
    norm = Normalize(vmin=min(y_true), vmax=max(y_true))
    df_pred['Color'] = df_pred['True'].apply(lambda x: cmap(norm(x)))

    plt.figure(figsize=(8, 6))
    sns.regplot(x='True', y='Pred', data=df_pred, scatter=False, color='black')
    sns.scatterplot(x='True', y='Pred', data=df_pred, palette=df_pred['Color'], hue='True', legend=False)

    plt.xlabel('True HAMD score')
    plt.ylabel('Predicted HAMD score')
    plt.title('ElasticNet Regression Correlation')

    plt.text(0.95, 0.25, f'R²: {r2:.3f}', ha='right', transform=plt.gca().transAxes)
    plt.text(0.95, 0.20, f'MAE: {mae:.3f}', ha='right', transform=plt.gca().transAxes)
    plt.text(0.95, 0.15, f'Pearson: {pearson_corr:.3f}', ha='right', transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def main():
    df = load_data('./mirna.csv')
    X, y = df.drop(columns=['HAMD']), df['HAMD']

    ratios = np.arange(0.01, 1.01, 0.05)
    alphas = np.concatenate(([10**-6, 10**-5, 10**-4, 10**-3, 10**-2], np.arange(0.1, 1.0, 0.1), np.arange(1, 11, 1)))

    best_ratio, best_alpha = elasticnet_grid_search(X, y, ratios, alphas)
    print(f"\nBest alpha: {best_alpha}, Best l1_ratio: {best_ratio}\n")

    feature_freq_df = bootstrap_feature_selection(X, y, best_alpha, best_ratio)
    feature_freq_df.to_csv('feature_selection_frequency.csv', index=False)

    best_thresh = find_best_threshold(X, y, feature_freq_df, best_alpha, best_ratio)
    print(f"\nBest threshold: {best_thresh}")

    mse, r2, mae, pearson_corr, y_true, y_pred = final_evaluation(X, y, best_alpha, best_ratio, best_thresh, feature_freq_df)

    print(f"\nFinal Evaluation:\nR²: {r2:.3f}, MAE: {mae:.3f}, Pearson corr: {pearson_corr:.3f}")

    plot_correlation(y_true, y_pred, r2, mae, pearson_corr, 'mirna_correlation_distribution.pdf')

if __name__ == '__main__':
    main()