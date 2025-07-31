import pandas as pd

def load_and_split_data(filepath):
    data_row = pd.read_csv(filepath, encoding='GBK')
    train_data = data_row[data_row.iloc[:, 0].str.startswith('train')]
    test_data = data_row[data_row.iloc[:, 0].str.startswith('test')]
    return train_data, test_data

def process_data(data, group_filter, replace_map, drop_columns):
    filtered_data = data[data['group'].isin(group_filter)].copy()
    filtered_data.loc[:, 'group'] = filtered_data['group'].replace(replace_map)
    group = filtered_data['group']
    features = filtered_data.drop(columns=drop_columns)
    return features, group

def prepare_datasets(filepath):
    train_data, test_data = load_and_split_data(filepath)

    drop_columns = ['allRNA', 'Hospital', 'Sample_id', 'Company', 'Batch', 
                    'group', 'Age', 'HAMD', 'Diagnosis', 'Gender']

    # MDD
    train_mdd_feature, train_mdd_group = process_data(train_data, [3, 1, 0], {3: 0}, drop_columns)
    test_mdd_feature, test_mdd_group = process_data(test_data, [3, 1, 0], {3: 0}, drop_columns)

    # BD
    train_bd_feature, train_bd_group = process_data(train_data, [2, 1], {2: 0}, drop_columns)
    test_bd_feature, test_bd_group = process_data(test_data, [2, 1], {2: 0}, drop_columns)

    return {
        'MDD': (train_mdd_feature, train_mdd_group, test_mdd_feature, test_mdd_group),
        'BD': (train_bd_feature, train_bd_group, test_bd_feature, test_bd_group)
}
