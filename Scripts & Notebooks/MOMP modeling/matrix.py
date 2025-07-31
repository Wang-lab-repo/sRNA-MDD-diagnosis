import seaborn as sns
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import numpy as np

for name, model in best_models.items():
    y_test_pred = model.predict(X_test)

    cm = confusion_matrix(y_test, y_test_pred)

    cm_percentage = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

    plt.figure(figsize=(8, 6))

    sns.heatmap(cm_percentage, annot=False, fmt=".2f", cmap='Blues', cbar=True, square=True,
                xticklabels=['Negative', 'Positive'], 
                yticklabels=['Negative', 'Positive'])

    for i in range(cm_percentage.shape[0]):
        for j in range(cm_percentage.shape[1]):
            color = 'black' if cm_percentage[i, j] < 50 else 'white'
            plt.text(j + 0.5, i + 0.5, f"{cm_percentage[i, j]:.2f}%", 
                     ha='center', va='center', color=color, fontsize=12)

    plt.title(f'{name} Confusion Matrix (Percentage)', fontsize=14)
    plt.xlabel('Predicted Label', fontsize=12)
    plt.ylabel('True Label', fontsize=12)

    plt.tight_layout()

    plt.savefig('{name}_confusion_matrix.svg', format='svg')

    plt.show()