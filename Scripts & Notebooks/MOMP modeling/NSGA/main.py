if __name__ == "__main__":
    data_path = './all.csv'
    dataset_parts = prepare_datasets(data_path)

    test_dataset = {
        "MDD": (dataset_parts['MDD'][2], dataset_parts['MDD'][3]),
        "BD": (dataset_parts['BD'][2], dataset_parts['BD'][3])
    }
    
    model_names = ['AdaBoost', 'CatBoost', 'GBDT', 'LightGBM', 
                   'Logistic Regression', 'MLP', 'Random Forest', 'SVM', 'XGBoost']

    models_dict = {
        "MDD": load_models(model_names, 'MDD'),
        "BD": load_models(model_names, 'BD')
    }

    performance_df = evaluate_model_performance(models_dict, test_dataset)
    performance_df.to_csv('./result/model_performance.csv', index=False)

    pareto_pool_df = optimize_model_selection(performance_df, ["AUC", "AUPRC", "Log-Loss", "InferenceTime"])