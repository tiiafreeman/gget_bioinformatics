import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

input_folder = "my_modules/analyzing_motif_matrix/input_folder"
output_folder = "my_modules/analyzing_motif_matrix/output_folder"
os.makedirs(output_folder, exist_ok=True)

matrix_file = os.path.join(input_folder, "VEEV_nsP3_motif_matrix.csv")
matrix_df = pd.read_csv(matrix_file, index_col=0)

motif_percentages = (matrix_df.sum() / len(matrix_df)) * 100
stats_df = pd.DataFrame({
    'Motif': motif_percentages.index,
    'Percentage': motif_percentages.values
})

stats_df['Prefix'] = stats_df['Motif'].str.split('_').str[0]

for prefix, group in stats_df.groupby('Prefix'):
    plt.figure(figsize=(10,6))
    group = group.sort_values('Percentage', ascending=False)
    sns.barplot(data=group, x='Motif', y='Percentage', palette='viridis')

    plt.title(f'Percentage of Isolates with {prefix} Motifs')
    plt.xlabel('Motif_ID')
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0,105)

    for i, v in enumerate(group['Percentage']):
        plt.text(i, v + 1, f'{v:.1f}%', ha='center', fontsize=9)

    plt.tight_layout()

    save_path = os.path.join(output_folder, f'{prefix}_motif_conservation.png')
    plt.savefig(save_path)
    print(f'Saved plot for {prefix} motifs to {save_path}')
    plt.close()