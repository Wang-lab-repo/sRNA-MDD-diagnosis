import joblib

def load_models(model_names, condition_tag):
    model_dict = {}
    for name in model_names:
        model_path = f'./model/{name}_model_{condition_tag}.pkl'
        try:
            model = joblib.load(model_path)
            model_dict[name] = model
            print(f"{name} model ({condition_tag}) loaded successfully.")
        except Exception as e:
            print(f"Failed to load {name} model ({condition_tag}): {e}")
    return model_dict
