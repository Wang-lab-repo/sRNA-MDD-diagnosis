import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'

plt.figure(figsize=(10, 8))
for name, (fpr, tpr, auc) in roc_curves.items():
    lower_bound, upper_bound = conf_intervals[name]
    plt.plot(fpr, tpr, label=f'{name} ( AUC = {auc:.3f} [{lower_bound:.3f}-{upper_bound:.3f}] )')

plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.01])
plt.xlabel('False Positive Rate', fontsize=12)
plt.ylabel('True Positive Rate', fontsize=12)
plt.title('ROC Curves', fontsize=15)
plt.legend(loc='lower right')

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.grid(False)

plt.savefig('ROC.svg', format='svg')

plt.show()