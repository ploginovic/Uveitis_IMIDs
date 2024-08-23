import datetime
import os
from lifelines import CoxPHFitter

import pandas as pd

from survival_analysis.backward_elimination import * 
from datetime import datetime

import pickle

def fit_and_show_cph(data, disease, grs, uve_col='uve_any',
                     age_col='age_uve_years', add_diseases=None, 
                     add_columns=['Sex_Female'], duration_col=None,
                     event_col=None, save_summary=True,
                     save_dir = None, return_model=True,
                     get_dummies = None):
    """
    Fit and display a summary of a Cox Proportional Hazards model using the provided data.
    Optionally, saves the summary to an Excel file.

    Parameters:
    data (pd.DataFrame): The dataset containing the patient data.
    disease (str): The base name of the disease column, appended with '_any'.
    grs (str): Column name for the genetic risk score.
    uve_col (str, optional): Column name indicating uveitis status. Defaults to 'uve_any'.
    age_col (str, optional): Column name for age at uveitis diagnosis. Defaults to 'age_uve_years'.
    add_columns (list, optional): Additional patient attributes to include. Defaults to ['Sex_Female'].
    duration_col (str, optional): Duration from uveitis diagnosis to the event. Defaults if None.
    event_col (str, optional): Event indicator column name. Defaults if None.
    save_summary (bool, optional): Flag to save the summary to an Excel file.

    Returns:
    None: The function prints and optionally saves the fit summary of the Cox model.
    """
    
    if add_columns is None:
        add_columns = []

    # Strip '_any' from the disease name to handle its normal form
    disease_name = disease.strip('_any')
    print(f"Disease formatted as {disease_name}")

    # Setup column names for the event and duration if not provided
    event_col = f'first_uve_{disease_name}' if event_col is None else event_col
    
    print(f"Event column is {event_col}")
    duration_col = f'uve_to_{disease_name}_years' if duration_col is None else duration_col
    print(f"Duration column is {duration_col}")

    # Define columns to export from the dataset
    cols_to_export = [grs, event_col, duration_col, age_col] + add_columns

    # Prepare data using the helper function
    cph_data = prepare_cph_data(data, disease, grs, uve_col,
                                age_col, add_columns, duration_col, event_col)
    

    if get_dummies is not None:
    # Note age groups that have no observations 
        gdf = cph_data.groupby(by=event_col)[get_dummies].value_counts().reset_index(name='counts')
        vals = gdf[(gdf[event_col]==1) & (gdf['counts']==0)][get_dummies].values

        cph_data = pd.get_dummies(cph_data, columns=[get_dummies], drop_first = False)
        # Selecting column with highest obs and dropping it
        most_occ = (cph_data[[i for i in cph_data.columns if get_dummies in i]]
                       .sum().idxmax())
        print(most_occ)
        cph_data.drop(most_occ, axis=1, inplace= True)

        # Dropping dummies with no observations
        cph_data.drop(columns=[i for i in cph_data.columns
                               if i.strip(f"{get_dummies}_") in vals.astype(str)],
                      inplace=True)
                  
    # Fit the Cox Proportional Hazards model and print the summary
    cph = CoxPHFitter()
    cph.fit(cph_data, duration_col=duration_col, event_col=event_col)
    cph.print_summary()
    
    
    if save_summary:
        # formatting filename of the saved xlsx file
        age_type = 'cont' if age_col=='age_uve_years' else age_col
        num_vars = len(cols_to_export)-2
        
        # Semiformatted file name
        semi_name = f"{disease_name}_{num_vars}_{age_type}"
        if save_dir is not None:
            save_cph_summary(cph, save_name=semi_name, directory=save_dir)
        else:
            save_cph_summary(cph, save_name=semi_name)
    
    if return_model:
        return cph, cph_data

        
def prepare_cph_data(data, disease, grs, uve_col='uve_any',
                     age_col='age_uve_years', add_columns=None, 
                     duration_col=None, event_col=None):
    """
    Prepare data for fitting a Cox Proportional Hazards model by setting up column names,
    filtering data, and selecting relevant columns.

    Parameters:
    data (pd.DataFrame): The dataset containing the patient data.
    disease (str): The base name of the disease column, appended with '_any'.
    grs (str): Column name for the genetic risk score.
    uve_col (str, optional): Column name indicating uveitis status. Defaults to 'uve_any'.
    age_col (str, optional): Column name for age at uveitis diagnosis in years. Defaults to 'age_uve_years'.
    add_columns (list, optional): Additional patient attributes to include in the analysis.
    duration_col (str, optional): Column name for the duration from uveitis diagnosis to the event. If None, a default based on disease is used.
    event_col (str, optional): Column name for the event indicator. If None, a default based on disease is used.

    Returns:
    pd.DataFrame: The filtered data ready for fitting a Cox model.
    """
    if add_columns is None:
        add_columns = []

    # Strip '_any' from the disease name to handle its normal form
    disease_name = disease.strip('_any')
    print(f"Disease formatted as {disease_name}")

    # Setup column names for the event and duration if not provided
    event_col = f'first_uve_{disease_name}' if event_col is None else event_col
    
    print(event_col)
    duration_col = f'uve_to_{disease_name}_years' if duration_col is None else duration_col
    
    print(duration_col)

    # Define columns to export from the dataset
    cols_to_export = [grs, event_col, duration_col, age_col] + add_columns

    # Filter the data for relevant cases
    cph_data = data.loc[(data[uve_col] == 1) & (data[f"first_{disease_name}"] != 1), cols_to_export]

    return cph_data

def save_cph_summary(cph_model, save_name,
                     directory='cph_summaries'):
    """
    Saves the summary of a CPH model to an Excel file with.

    Parameters:
    cph_model (CoxPHFitter object): Fitted Cox Proportional Hazards model.
    disease_name (str): Name of the disease used in the filename.

    Returns:
    None
    """
    
    # Create the directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Format today's date as YYYYMMDD
    today = datetime.date.today().strftime('%d%m%Y')

    # Create filename with today's date and the disease name
    filename = f"{save_name}_CPH_summary_{today}.xlsx"
    
    file_path = os.path.join(directory, filename)

    # Save the model summary to an Excel file
    cph_model.summary.round(4).to_excel(file_path)
    print(f"Saved Cox PH model summary to {file_path}")
    
    
def prepare_and_save_model(data, disease, grs, event_col, duration_col,
                           age_col, add_columns, save_path, return_model=True,
                           return_data=True,
                           trace=False, penalizer = 0):
    """
    Prepare data, perform backward elimination, check model assumptions, and save the model with a date and disease name in the filename.

    Parameters:
        data (pd.DataFrame): Input dataset.
        disease (str): Disease name for the data preparation and filename.
        grs (str): Genetic Risk Score column name.
        event_col (str): Event column name.
        duration_col (str): Duration column name.
        age_col (str): Age column name.
        add_columns (list): List of additional columns to include.
        save_path (str): Path to save the finalized model.

    Returns:
        None: Saves the model to the specified path and prints the path.
    """
    # Prepare the data
    preped_df = prepare_cph_data(data, disease=disease, grs=grs, event_col=event_col,
                                 duration_col=duration_col, age_col=age_col, add_columns=add_columns)

    # Perform backward elimination based on p-values
    features, cph = backward_elimination_p(preped_df,
                                           duration_col=duration_col,
                                           event_col=event_col,
                                           penalizer=penalizer)
    
    print(f'Backward elimination based on P-value completed.'
          +f'Features reminaing: {features} \n')
    
    
    # Perform backward elimination based on partial AIC
    
    features, cph = backward_elimination_AIC(preped_df[features 
                                                       + [duration_col, event_col]],
                                             duration_col=duration_col,
                                             event_col=event_col)

    # Check model assumptions
    assumption_check = cph.check_assumptions(preped_df[features 
                                                       + [duration_col, event_col]],
                                             p_value_threshold=0.05)
    if trace:
        print("Assumption Check Results:", assumption_check, '\n')

    # Format the current date
    current_date = datetime.now().strftime("%Y-%m-%d")

    # Save the model
    filename = f'{save_path}/{disease}_finalised_model_{current_date}.pkl'
    with open(filename, 'wb') as file:
        pickle.dump(cph, file)
    print(f"Model saved to {filename}")
    
    if return_model and not return_data:
        return cph
    elif return_model and return_data:
        return cph, preped_df[features + [duration_col, event_col]]
    
    elif return_data:
        preped_df[features + [duration_col, event_col]]
    