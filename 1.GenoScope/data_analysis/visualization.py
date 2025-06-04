import matplotlib.pyplot as plt
import seaborn as sns

def plot_correlation_matrix(df, output_path=None):
    """
    Строит и (опционально) сохраняет корреляционную матрицу.
    """
    try:
        corr = df.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap='coolwarm')
        plt.title("Correlation Matrix")
        if output_path:
            plt.savefig(output_path)
        plt.show()
    except Exception as e:
        print(f"Error plotting correlation matrix: {e}")

def plot_distributions(df, output_dir=None):
    """
    Строит гистограммы распределения для всех числовых признаков.
    """
    try:
        numeric_cols = df.select_dtypes(include=['float64', 'int']).columns
        for col in numeric_cols:
            plt.figure()
            sns.histplot(df[col], kde=True)
            plt.title(f"Distribution of {col}")
            if output_dir:
                plt.savefig(f"{output_dir}/{col}_distribution.png")
            plt.show()
    except Exception as e:
        print(f"Error plotting distributions: {e}")

def plot_pca(pca_result, target_column=None):
    """
    Визуализирует результаты PCA.
    """
    try:
        plt.scatter(pca_result['PC1'], pca_result['PC2'], c=target_column, cmap='viridis')
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA Visualization')
        if target_column is not None:
            plt.colorbar(label="Target")
        plt.show()
    except Exception as e:
        print(f"Error plotting PCA: {e}")
