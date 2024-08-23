from lifelines import CoxPHFitter
import numpy as np
import pandas as pd
from itertools import combinations

def create_dummies(data):
    for col in data:
        print(col)
        if data[col].isnull().values.any() or data[col].dtype =='O':
            print(col, ' needs a dummy variable: empty data')

            column_to_drop = str(col)+"_"+str(data[col].value_counts().idxmax())
            print("The column will be dropped: ",column_to_drop)
            data = pd.get_dummies(data, columns = [col], drop_first=False)
            data = data.loc[:, data.columns!=str(column_to_drop)]
            
            
        if data[col].dtype=='int64' and (data[col].min!=0 or data[col].max!=1):
            print(col, ' needs a dummy variable: categorical intigers suspected')
            
            column_to_drop = str(col)+"_"+str(0)
            print("The column will be dropped: ",column_to_drop)
            data = pd.get_dummies(data, columns = [col], drop_first=False)
            data = data.loc[:, data.columns!=str(column_to_drop)]  
            
    return(data)


def backward_elimination_AIC(data, duration_col='uve_to_MS_years', event_col='first_uve_MS'):
    """
    Perform backward elimination by first checking the model with all features and then testing all possible sub-combinations 
    of features at each step based on AIC, always retaining features that contain 'age' and 'sex' in their names.

    Parameters:
        data (DataFrame): The dataset containing all necessary columns.
        duration_col (str): The name of the column representing the duration until the event.
        event_col (str): The name of the column representing the event occurrence (binary).

    Returns:
        list, CoxPHFitter: The best combination of features and the corresponding fitted model.
    """
    # Identify features that must always be retained (features containing 'age' or 'sex')
    mandatory_features = [col for col in data.columns if 'age' in col.lower() or 'sex' in col.lower()]
    features = [col for col in data.columns if col not in [event_col, duration_col] + mandatory_features]

    # Start with all features, including mandatory ones
    all_features = features + mandatory_features

    # Fit the initial model with all candidate features
    initial_model = CoxPHFitter()
    initial_model.fit(data[all_features + [duration_col, event_col]], duration_col=duration_col, event_col=event_col)
    best_AIC = initial_model.AIC_partial_
    best_features = all_features
    best_model = initial_model

    print(f"Initial model with all features AIC: {best_AIC}")

    # Iteratively remove one feature at a time
    while features:
        current_best_AIC = float('inf')
        current_best_combo = None
        current_best_model = None

        # Generate all combinations of features with one less feature at each step
        for combo in combinations(features, len(features) - 1):
            current_features = list(combo) + mandatory_features + [duration_col, event_col]
            model = CoxPHFitter()
            model.fit(data[current_features], duration_col=duration_col, event_col=event_col)
            current_AIC = model.AIC_partial_

            # Find the best combination at this level
            if current_AIC < current_best_AIC:
                current_best_AIC = current_AIC
                current_best_combo = combo
                current_best_model = model

        # Compare the best combination of the current step with the overall best
        if current_best_AIC + 2 < best_AIC:
            best_AIC = current_best_AIC
            features = list(current_best_combo)  # Reduce the feature set for the next iteration
            best_features = features + mandatory_features
            best_model = current_best_model
            print(f"Reduced to features {best_features} with AIC: {best_AIC}")
        else:
            print(f"No further improvement;"
                  +f" best acheived was {current_best_AIC}\n"+
                  "stopping elimination.")
            break

    return best_features, best_model


#    return best_features, best_model

def backward_elimination_p(data, duration_col='uve_to_MS_years', event_col='first_uve_MS', significance_level=0.05,
                          penalizer = 0.002):
    """
    Perform backward elimination based on p-values from a Cox Proportional Hazards model,
    ensuring features related to 'age' and 'sex' are always retained.

    Parameters:
        data (pd.DataFrame): The dataset containing the patient data.
        duration_col (str): Column name for the duration from event diagnosis to the event.
        event_col (str): Column name for the event indicator (1 if event occurred, 0 otherwise).
        significance_level (float): The p-value threshold above which a feature should be considered for removal.

    Returns:
        tuple: A tuple containing the list of remaining features after elimination and the final Cox model.
    """
    # Identify features that must always be retained (features containing 'age' or 'sex')
    mandatory_features = [col for col in data.columns if 'age' in col.lower() or 'sex' in col.lower()]

    # Initial set of features excludes the duration and event columns and mandatory features
    features = [col for col in data.columns if col not in [event_col, duration_col] + mandatory_features]
    all_features = mandatory_features + features  # Full feature set including mandatory features

    while len(features) > 0:
        model = CoxPHFitter(penalizer=penalizer).fit(data[all_features + [duration_col, event_col]], duration_col=duration_col, event_col=event_col)
        p_values = model.summary['p']

        # Exclude mandatory features from the removal candidates
        removable_p_values = p_values.drop(index=mandatory_features, errors='ignore')

        max_p_value = removable_p_values.max()
        if max_p_value >= significance_level:
            worst_feature = removable_p_values.idxmax()  # Get the feature with the maximum p-value
            features.remove(worst_feature)  # Remove the feature with the highest p-value from features
            all_features.remove(worst_feature)  # Also remove it from all_features
            print(f"Removed {worst_feature} with p-value: {max_p_value}")
        else:
            print("No further improvements; p-values are below the significance level.")
            break
    
    # Fit the final model with the retained features
    final_model = CoxPHFitter().fit(data[all_features + [duration_col, event_col]], duration_col=duration_col, event_col=event_col)
    return all_features, final_model



def stepwise_selection(data, target,SL_in=0.05,SL_out = 0.05):
    initial_features = data.columns.tolist()
    best_features = []
    while (len(initial_features)>0):
        remaining_features = list(set(initial_features)-set(best_features))
        new_pval = pd.Series(index=remaining_features)
        for new_column in remaining_features:
            model = sm.OLS(target, sm.add_constant(data[best_features+[new_column]])).fit()
            new_pval[new_column] = model.pvalues[new_column]
        min_p_value = new_pval.min()
        if(min_p_value<SL_in):
            best_features.append(new_pval.idxmin())
            while(len(best_features)>0):
                best_features_with_constant = sm.add_constant(data[best_features])
                p_values = sm.OLS(target, best_features_with_constant).fit().pvalues[1:]
                max_p_value = p_values.max()
                if(max_p_value >= SL_out):
                    excluded_feature = p_values.idxmax()
                    best_features.remove(excluded_feature)
                else:
                    break 
        else:
            break
    return best_features


# def backward_elimination_AIC(data, significance_level = 0.05):
    
#     data=create_dummies(data)
    
#     features = data.columns.tolist()
#     print(" Initial variable list in the model: ", data.loc[:,((data.columns!='first_ON') & (data.columns!='ON_to_MS_years'))].columns.tolist(), '\n')


#     while(len(features)>0):
        
#         AIC_dict = {}
#         model = CoxPHFitter().fit(data, duration_col='ON_to_MS_years', event_col='first_ON')
#         current_AIC = model.AIC_partial_
        
#         for i in features:
#             if (i!='first_ON') and i!='ON_to_MS_years':
            
#                 new_model = CoxPHFitter().fit(data.loc[:, data.columns !=i], duration_col='ON_to_MS_years', event_col='first_ON')
#                 AIC_dict[i] = round(new_model.AIC_partial_,3)               

#         print('AIC after removing each variable:', AIC_dict)    
#         min_new_AIC= min(AIC_dict.values())
#         removed_var = list(AIC_dict.keys())[list(AIC_dict.values()).index(min_new_AIC)]
        
#         if min_new_AIC < current_AIC: 
            
#             print('Current model: ')
#             #model.print_summary()
            
#             features.remove(str(removed_var))            
#             current_AIC = min_new_AIC
#             print('\n','Removed variable: ', removed_var, "\n", "New AIC is ", min_new_AIC, '\n')
#             data=data.loc[:, features]
            
#         elif min_new_AIC >=current_AIC:
#             print("\n" , "Model optimised: no lower AIC can be achieved")
#             print("The final set of variables is: ", features)
#             model = CoxPHFitter().fit(data, duration_col='ON_to_MS_years', event_col='first_ON')
#             print("\n")
#             #model.print_summary()
#             break

#     return features, model
# def backward_elimination_p(data, significance_level = 0.05):
#     features = data.columns.tolist()
    
#     while(len(features)>0):
#         model = CoxPHFitter().fit(data, duration_col='ON_to_MS_years', event_col='first_ON')
#         p_values = model.summary.loc[:,'p']
#         max_p_value = p_values.max()
        
#         if(max_p_value >= significance_level):
#             excluded_feature = CoxPHFitter().summary[(CoxPHFitter().summary.p == max(CoxPHFitter().summary.p))].index[0]
#             features.remove(str(excluded_feature))
#             data = data.loc[:,features]
#         else:
#             break 
#     return features