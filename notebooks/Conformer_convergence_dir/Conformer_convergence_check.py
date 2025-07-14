#!/usr/bin/env python
# coding: utf-8

# # check the converge of the HL-gap wrt. to the number of conformers

# In[1]:


import sys, os
import pandas as pd


# In[2]:


notebook_dir = os.path.abspath('')  # current directory: /home/sat/HL-gaps-pub/notebooks
project_root = os.path.abspath(os.path.join(notebook_dir, '..'))  # /home/sat/HL-gaps-pub
sys.path.insert(0, project_root)


# In[3]:


from hl_gaps_pub.hl_gaps_pub import calculate_gap, embed_confs


# In[4]:


try:
    df = pd.read_pickle("../../data/processed/gap_smile.pkl")
except FileNotFoundError:
    print("Error: 'data/processed/gap_smile.pkl' not found.")


# In[5]:


#df


# In[6]:


# The functions calculate_gap and embed_confs are assumed to be imported
# from hl_gaps_pub.hl_gaps_pub as shown in your notebook.

def analyze_conformer_convergence(df, smiles_column='SMILE', num_smiles_to_sample=10, conformer_counts=[5, 10, 20, 40]):
    """
    Calculates HL-gaps for a random selection of SMILES from a DataFrame,
    varying the number of conformers.

    Args:
        df (pd.DataFrame): DataFrame containing SMILES strings.
        smiles_column (str): Name of the column in df with SMILES strings.
                             Defaults to 'SMILES'.
        num_smiles_to_sample (int): Number of random SMILES to select from df.
                                    Defaults to 10.
        conformer_counts (list): List of integers representing the number of
                                 conformers to generate for each SMILES.
                                 Defaults to [5, 10, 20, 50].

    Returns:
        list: A list of dictionaries, where each dictionary contains:
              'smiles': The SMILES string.
              'num_conformers': The number of conformers used.
              'hl_gap': The calculated HL-gap (float or None if error).
              'error': Error message if calculation failed (str or None).
    """
    results = []

    if not isinstance(df, pd.DataFrame):
        print("Error: Input 'df' is not a pandas DataFrame.")
        return results
    if smiles_column not in df.columns:
        print(f"Error: Smiles column '{smiles_column}' not found in DataFrame.")
        return results
    if df.empty:
        print("Error: Input DataFrame 'df' is empty.")
        return results

    # Ensure num_smiles_to_sample is not greater than available SMILES
    actual_smiles_to_sample = min(num_smiles_to_sample, len(df[smiles_column].unique()))
    
    if actual_smiles_to_sample == 0:
        print("No SMILES available to sample.")
        return results

    # Sample unique SMILES to avoid redundant processing if duplicates exist
    # and sampling is less than total unique SMILES.
    # If num_smiles_to_sample is larger than unique SMILES, it will take all unique ones.
    sampled_smiles = df[smiles_column].drop_duplicates().sample(
        n=actual_smiles_to_sample,
        random_state=42  # for reproducibility
    ).tolist()

    print(f"Selected {len(sampled_smiles)} SMILES for analysis: {sampled_smiles}")

    for smile_string in sampled_smiles:
        print(f"\nProcessing SMILES: {smile_string}")
        for n_confs in conformer_counts:
            print(f"  Attempting with {n_confs} conformers...")
            hl_gap_value = None
            error_message = None
            try:
                # Step 1: Embed conformers
                # embed_confs is expected to return a molecule object with conformers,
                # or None/raise an exception on failure.
                mol_with_confs = embed_confs(smile=smile_string, num_confs=n_confs)

                if mol_with_confs is None:
                    error_message = f"Conformer generation failed for {n_confs} conformers."
                    print(f"    {error_message}")
                else:
                    # Step 2: Calculate HL-gap
                    # calculate_gap is expected to take the output of embed_confs
                    hl_gap_value = calculate_gap(molecule=mol_with_confs, method='GFN2-xTB', accuracy=1.0, temperature=300.0)
                    if hl_gap_value is not None: # Assuming calculate_gap might return None on error
                        print(f"    Successfully calculated HL-gap: {hl_gap_value:.4f} eV (for {n_confs} conformers)")
                    else:
                        error_message = "HL-gap calculation returned None."
                        print(f"    {error_message}")


            except Exception as e:
                error_message = f"An error occurred: {str(e)}"
                print(f"    {error_message}")
            
            results.append({
                'smiles': smile_string,
                'num_conformers': n_confs,
                'hl_gap': hl_gap_value,
                'error': error_message
            })
            
    return results


# In[ ]:


results= analyze_conformer_convergence(df)


# In[ ]:


print(results)
#results.save


# In[ ]:


for result in results:
    print(result)

# Or convert results to a DataFrame for easier analysis:
results_df = pd.DataFrame(results)
#results_df.head()


# In[ ]:


# --- Check if results_df exists and is not empty before saving ---
if 'results_df' in locals() and isinstance(results_df, pd.DataFrame) and not results_df.empty:
    # 1. Save to CSV
    csv_file_path = "conformer_convergence_results.csv"
    results_df.to_csv(csv_file_path, index=False)  # index=False prevents pandas from writing the DataFrame index as a column
    print(f"Results saved to CSV: {csv_file_path}")

    # 2. Save to Pickle
    pickle_file_path = "conformer_convergence_results.pkl"
    results_df.to_pickle(pickle_file_path)
    print(f"Results saved to Pickle: {pickle_file_path}")
else:
    print("Error: 'results_df' is not defined, not a DataFrame, or is empty. Cannot save.")

# --- To load them back later (for verification or future use) ---
# Load from CSV
# loaded_df_csv = pd.read_csv("conformer_convergence_results.csv")
# print("\\nLoaded from CSV:")
# print(loaded_df_csv.head())

# Load from Pickle
# loaded_df_pickle = pd.read_pickle("conformer_convergence_results.pkl")
# print("\\nLoaded from Pickle:")
# print(loaded_df_pickle.head())

