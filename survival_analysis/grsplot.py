import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

from statannot import add_stat_annotation


def define_groups(data, disease_col, uveitis_col, disease_label):
    """
    Defines groups based on the presence of a disease and uveitis.

    Parameters:
    data (pd.DataFrame): The DataFrame containing patient data.
    disease_col (str): The column name for binary disease.
    uveitis_col (str): The column name for binary uveitis.
    disease_label (str): A human-readable label for the disease.

    Returns:
    dict: A dictionary containing DataFrames for each group.
    """
    # Controls are cases without the disease and without uveitis
    controls = data.loc[(data[disease_col] != 1) & (data[uveitis_col] != 1)]

    # Disease only are cases with the disease but without uveitis
    disease_only = data.loc[(data[disease_col] == 1) & (data[uveitis_col] != 1)]

    # Uveitis only are cases without the disease but with uveitis
    uve_only = data.loc[(data[disease_col] != 1) & (data[uveitis_col] == 1)]

    # Disease-Uve are cases with both the disease and uveitis
    disease_uve = data.loc[(data[disease_col] == 1) & (data[uveitis_col] == 1)]

    return {
        'Controls': controls,
        f'{disease_label} only': disease_only,
        'Uve': uve_only,
        f'{disease_label}-Uve': disease_uve
    }


def plot_distributions(ax, groups, grs_col, colors, disease_name):
    """
    Plots Kernel Density Estimation (KDE) for each group on the given axis using specified colors.

    Parameters:
    ax (matplotlib.axes.Axes): The axes object to plot the distributions.
    groups (dict): Dictionary of groups with their respective data.
    grs_col (str): Column name of the genetic risk scores.
    colors (dict): Dictionary mapping group names to colors for the plot.
    disease_name (str): The name of the disease derived from the column name.
    """
    for group_name, group_data in groups.items():
        
        label = f"{group_name}, n={len(group_data):,}"
        sns.kdeplot(group_data[grs_col], ax=ax, color=colors[group_name], label=label)


def perform_tests(groups, grs_col, comparisons):
    p_values = []
    for group1, group2 in comparisons:
        result = stats.ttest_ind(groups[group1][grs_col],
                                 groups[group2][grs_col],
                                 equal_var=False)
        p_values.append(result.pvalue)
    return p_values

def correct_p_values(p_values, comparisons):
    return [p for p in p_values]

def annotate_results(ax, comparisons, p_values):
    """
    Adds annotations to a plot with the results of statistical tests, placing them in the bottom left with a 
    semitransparent background for better readability. Each p-value is formatted to show two significant digits or special notation for very small values.

    Parameters:
    ax (matplotlib.axes.Axes): The axes object where annotations should be added.
    comparisons (list of tuple): List of group comparisons.
    p_values (list): List of corrected p-values associated with the comparisons.
    """
    text_str = ""
    for (group1, group2), p in zip(comparisons, p_values):
        formatted_p = format_p_value(p)  # Use the helper function to format the p-value appropriately
        text_str += f'{group1} vs {group2}: P{formatted_p}\n'
    
    # Using bbox to create a semitransparent white background for the text
    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="b", lw=0.2, alpha=0.8)
    # Position the text at the bottom left of the axes (transform=ax.transAxes to use axis coordinates)
    ax.text(0.05, 0.05, text_str.strip(), transform=ax.transAxes, verticalalignment='bottom', 
            horizontalalignment='left', fontsize=8, bbox=bbox_props)


def format_p_value(p):
    """
    Formats the p-value to display with two significant digits or returns '<0.0001' for p-values smaller than 0.0001.

    Parameters:
    p (float): The p-value to format.

    Returns:
    str: Formatted p-value, ensuring two significant digits are displayed, or a special notation for very small values.
    """
    if p > 1:
        return "=1"
    if p < 0.0001:
        return "<0.0001"
    elif p < 0.001 :
        return f"={p:.1g}"  # Use scientific notation for small numbers not covered by the first case
    else:
        return f"={p:.2g}"  # Use general format for numbers larger than or equal to 0.001



def plot_grs(data, disease_col, grs_col, uveitis_col, ax, colors=None,
             fontsize=10):
    """
    Main function to plot genetic risk scores and perform statistical analysis.
    Dynamically determines disease names based on the provided disease_col by stripping '_any'.

    Parameters:
    data (pd.DataFrame): The dataset containing the data.
    disease_col (str): Column name for the disease status.
    grs_col (str): Column name for the genetic risk scores.
    uveitis_col (str): Column name for uveitis status.
    ax (matplotlib.axes.Axes): The axes object to plot the KDE and annotations.
    colors (dict, optional): Dictionary mapping group names to colors for the plot.
    """
    disease_name = disease_col.replace('_any', '')
    

    if colors is None:
        colors = {
            'Controls': 'red',
            f'{disease_name} only': 'blue',
            'Uve': 'green',
            f'{disease_name}-Uve': 'purple'
        }

    groups = define_groups(data, disease_col, uveitis_col, disease_name)
    plot_distributions(ax, groups, grs_col, colors, disease_name)
    ax.set_title(f'Distribution of {disease_name} GRS')
    ax.set_xlim(None, 8)
    ax.legend(fontsize=fontsize)

    comparisons = [
        ('Controls', f'Uve'),
        (f'{disease_name} only', 'Uve'),
        ('Uve', f'{disease_name}-Uve'),
        (f'{disease_name}-Uve', f'{disease_name} only')
    ]
    
    p_values = perform_tests(groups, grs_col, comparisons)
    # corrected_p_values = correct_p_values(p_values, comparisons)
    corrected_p_values = [float(i)*(8*4) for i in p_values]
    annotate_results(ax, comparisons, corrected_p_values)
    
    

def plot_violin(data, disease_col, grs_col, uveitis_col, ax, colors=None,
                alpha_values=None,fontsize=10):
    disease_name = disease_col.replace('_any', '')
    groups = define_groups(data, disease_col, uveitis_col, disease_name)
    group_data = (pd.concat(groups.values(), keys=groups.keys())
                  .reset_index(level=0)
                  .rename(columns={'level_0': 'Group'})
                 )
    group_data.rename(columns={group_data.columns[1]: grs_col}, inplace=True)
    
    print(group_data.Group.value_counts())

    if not colors:
        unique_groups = group_data['Group'].unique()
        palette = sns.color_palette("hsv", len(unique_groups))
        colors = dict(zip(unique_groups, palette))

    order = ["Controls", "Uve",
             f"{disease_name}-Uve", f"{disease_name} only"]
    sns.violinplot(x='Group', y=grs_col, data=group_data,
                  inner_kws=dict(box_width=15, whis_width=2, color=".8"), palette=colors, ax=ax,
                   order=order)
    ax.set_xlabel('Group', fontsize=fontsize)
    ax.set_ylabel(grs_col, fontsize=fontsize)
    ax.set_title(f"Violin Plot of {disease_name} GRS", fontsize=fontsize)

    # Check groups and print them for debugging
    print("Unique groups in data:", group_data['Group'].unique())

    if alpha_values:
        for violin, alpha in zip(ax.collections[::2], alpha_values):
            violin.set_alpha(alpha)

    # Adjust box_pairs according to your actual group names
    box_pairs = [
        (f"{disease_name}-Uve", f"{disease_name} only"),
        ('Uve', f'{disease_name}-Uve'),
        ('Uve', "Controls"),
        ('Uve', f"{disease_name} only")
    ]

    # Ensure box_pairs are valid
    valid_pairs = [pair for pair in box_pairs if all(item in group_data['Group'].unique() for item in pair)]

    if not valid_pairs:
        print("Warning: No valid box pairs found!")
    else:
        add_stat_annotation(ax, data=group_data, x='Group', y=grs_col, order=order,
                            box_pairs=box_pairs,
                            test='t-test_welch', pvalue_format_string="{.3f}",
                            loc='inside',
                            fontsize=fontsize, comparisons_correction=None )


def plot_grs_pair(data, grs1, grs2, disease1='MS_any', disease2='SLE_any',
                  a1=0.1, a2=0.4, uveitis=False, filter_col='uve_any',
                  additional_filter=None):
    """
    Plot Genetic Risk Scores (GRS) for different diseases with optional filtering, including KDE plots on the axes.
    Perform Welch's t-test between the groups.

    Parameters:
    data (pd.DataFrame): DataFrame containing the data.
    grs1 (str): Column name for the first genetic risk score.
    grs2 (str): Column name for the second genetic risk score.
    disease1 (str): Column name for the first disease indicator.
    disease2 (str): Column name for the second disease indicator.
    a1 (float): Alpha value for the 'other cases' points.
    a2 (float): Alpha value for the disease points.
    uveitis (bool): Whether to filter data for uveitis cases.
    filter_col (str): Column name for the uveitis filter.
    additional_filter (dict): Additional filter to apply, specified as a dictionary with column names as keys and filter values as values.
    """
    # Apply uveitis filter if specified
    if uveitis:
        data = data.loc[data[filter_col] == 1]
    
    # Apply additional filter if specified
    if additional_filter:
        for col, value in additional_filter.items():
            data = data.loc[data[col] == value]

    disease1_cases = data[data[disease1] == 1]
    disease2_cases = data[data[disease2] == 1]
    other_cases = data[(data[disease1] != 1) & (data[disease2] != 1)]

    # Create a JointGrid for the scatter plot with KDE plots on the axes
    g = sns.JointGrid(data=data, x=grs1, y=grs2, space=0)

    # Plot each group with different colors
    g.ax_joint.scatter(other_cases[grs1], other_cases[grs2], color='green', label=f'Other Cases (n={len(other_cases)})', alpha=a1)
    g.ax_joint.scatter(disease2_cases[grs1], disease2_cases[grs2], color='orange', label=f'{disease2.replace("_any", "")} (n={len(disease2_cases)})', alpha=a2)
    g.ax_joint.scatter(disease1_cases[grs1], disease1_cases[grs2], color='blue', label=f'{disease1.replace("_any", "")} (n={len(disease1_cases)})', alpha=a2)

    # Add KDE plots to the marginal axes
    sns.kdeplot(data=disease1_cases[grs1], ax=g.ax_marg_x, color='blue', fill=True, alpha=0.3)
    # sns.kdeplot(y=disease1_cases[grs2], ax=g.ax_marg_y, color='blue', fill=True, alpha=0.3)
    # sns.kdeplot(data=disease2_cases[grs1], ax=g.ax_marg_x, color='orange', fill=True, alpha=0.3)
    sns.kdeplot(y=disease2_cases[grs2], ax=g.ax_marg_y, color='orange', fill=True, alpha=0.3)
    sns.kdeplot(data=other_cases[grs1], ax=g.ax_marg_x, color='green', fill=True, alpha=0.3)
    sns.kdeplot(y=other_cases[grs2], ax=g.ax_marg_y, color='green', fill=True, alpha=0.3)

    # Adding titles and labels
    g.ax_joint.set_title('Scatter Plot of Genetic Risk Scores' + (' in Uveitis' if uveitis else ''))
    g.ax_joint.set_xlabel(grs1)
    g.ax_joint.set_ylabel(grs2)
    g.ax_joint.legend()
    g.ax_joint.grid(True)

    # Perform Welch's t-test
    ttest1 = stats.ttest_ind(disease1_cases[grs1], other_cases[grs1], equal_var=False)
    ttest2 = stats.ttest_ind(disease2_cases[grs1], other_cases[grs1], equal_var=False)
    
    print(f"Welch's t-test results for {disease1.replace('_any', '')} vs Other Cases on {grs1}:")
    print(f"t-statistic: {ttest1.statistic:.4f}, p-value: {ttest1.pvalue:.4f}")
    
    # print(f"Welch's t-test results for {disease2.replace('_any', '')} vs Other Cases on {grs1}:")
    # print(f"t-statistic: {ttest2.statistic:.4f}, p-value: {ttest2.pvalue:.4f}")

    
    ttest3 = stats.ttest_ind(disease1_cases[grs2], other_cases[grs2], equal_var=False)
    ttest2 = stats.ttest_ind(disease2_cases[grs2], other_cases[grs2], equal_var=False)
    
    # print(f"Welch's t-test results for {disease1.replace('_any', '')} vs Other Cases on {grs1}:")
    # print(f"t-statistic: {ttest1.statistic:.4f}, p-value: {ttest1.pvalue:.4f}")
    
    print(f"Welch's t-test results for {disease2.replace('_any', '')} vs Other Cases on {grs2}:")
    print(f"t-statistic: {ttest2.statistic:.4f}, p-value: {ttest2.pvalue:.4f}")
    plt.show()