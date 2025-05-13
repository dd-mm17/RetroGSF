import pandas as pd
import numpy as np
import os


def rank_similar_solvents(target_smiles, data_path='RetroGSF/data/Solvant_properties_with_smile.csv', n_recommendations=5):
    """
    Find solvents with similar physical properties to the target solvent and rank them.
    
    Parameters:
    -----------
    target_smiles : str
        SMILES string of the solvent to find alternatives for
    data_path : str
        Path to the CSV file containing solvent data
    n_recommendations : int
        Number of recommendations to return (default: 5)
        
    Returns:
    --------
    dict of pandas DataFrames containing ranked alternatives or error message if SMILES not found
    """
    # Load the data
    df = pd.read_csv(data_path, sep=';')
    
    # Check if SMILES exists in the dataset
    if target_smiles not in df['SMILES'].values:
        return f"Error: No solvent found with SMILES '{target_smiles}' in the dataset."
    
    # Get the solvent name and properties
    target_solvent = df[df['SMILES'] == target_smiles]['Name'].iloc[0]
    target = df[df['SMILES'] == target_smiles].iloc[0]
    
    # Check if target solvent is hazardous and warn user
    if target['Adjusted ranking'] in ['Hazardous', 'Highly Hazardous']:
        print(f"Warning: '{target_solvent}' (SMILES: {target_smiles}) is classified as {target['Adjusted ranking']}. "
              "Recommendations will be for safer alternatives only.")
    
    # Create working copy of DataFrame for filtering
    df_working = df.copy()
    
    # Convert ranking columns to numeric at the start to avoid any string comparisons
    ranking_columns = ['Environment Ranking', 'Health Ranking', 'Safety Ranking']
    for column in ranking_columns:
        df_working[column] = pd.to_numeric(df_working[column], errors='coerce')
    
    # Filtering the data
    #Hazardous
    df_filtered = df_working[~df_working['Adjusted ranking'].isin(['Hazardous', 'Highly Hazardous'])].copy()
    
    #Too high envirement/heath/safety ranking
    mask = (
        (df_filtered['Environment Ranking'] <= 5) &
        (df_filtered['Health Ranking'] <= 5) &
        (df_filtered['Safety Ranking'] <= 5)
    )
    df_filtered = df_filtered[mask].copy()
    
    if len(df_filtered) == 0:
        return "Error: No solvents found meeting safety and environmental criteria (scores ≤ 5)."
    
    # Physical state outside the target solvant
    target_melting = target['Melting point']
    target_boiling = target['Boiling point']
    
    df_filtered['temp_range_compatible'] = (
        (df_filtered['Melting point'] <= target_melting + 15) &  # Melting point not more than 15°C higher
        (df_filtered['Boiling point'] >= target_boiling - 15)    # Boiling point not more than 15°C lower
    )
    
    # Filter based on temperature criteria
    df_filtered = df_filtered[df_filtered['temp_range_compatible']].copy()
    
    # If no solvents match the temperature criteria, inform the user
    if len(df_filtered) == 0:
        return f"Warning: No safe solvents found within ±15°C range"
    
    # Calculate similarity scores based on physical properties
    # The properties tested and their weights

    properties = {
        'Density': 0.1,
        'Dielectric': 0.4,
        'Dipole': 0.35,  
        'Refractive Index': 0.15,  
    }
    
    # Initialize similarity score
    df_filtered['similarity_score'] = 0
    
    # Calculate normalized distances for each property
    for prop, weight in properties.items():
        if pd.notna(target[prop]) and target[prop] != 'NA':
            # Convert values to numeric, treating 'NA' as NaN
            df_filtered[prop] = pd.to_numeric(df_filtered[prop], errors='coerce')
            
            # Skip if there's no valid data for this property
            if df_filtered[prop].isna().all():
                continue
            
            # Calculate the range of the property for normalization
            prop_range = df_filtered[prop].max() - df_filtered[prop].min()
            
            # Avoid division by zero
            if prop_range == 0:
                continue
            
            # Calculate normalized distance
            df_filtered[f'{prop}_distance'] = np.abs(df_filtered[prop] - float(target[prop])) / prop_range
            
            # Add to similarity score (weighted)
            df_filtered['similarity_score'] += weight * df_filtered[f'{prop}_distance']
    
    # Filter out the target solvent itself from recommendations
    df_similar = df_filtered[df_filtered['Name'] != target_solvent].copy()
    
    # Sort by similarity score (lower is more similar)
    df_similar = df_similar.sort_values('similarity_score')
    
    # Get top N most similar solvents
    top_similar = df_similar.head(n_recommendations).copy()
    
    # Create rankings
    # For each ranking, lower numbers are better
    
    ranking_columns = ['Environment Ranking', 'Health Ranking', 'Safety Ranking']
    for column in ranking_columns:
        top_similar.loc[:, column] = pd.to_numeric(top_similar[column], errors='coerce')
    
    # Calculate overall ranking as a weighted sum (equal weights)
    top_similar['Overall Ranking'] = (
        top_similar['Environment Ranking'] + 
        top_similar['Health Ranking'] + 
        top_similar['Safety Ranking']
    ) / 3
    
    # Sort by various criteria
    env_ranked = top_similar.sort_values('Environment Ranking')
    health_ranked = top_similar.sort_values('Health Ranking')
    safety_ranked = top_similar.sort_values('Safety Ranking')
    overall_ranked = top_similar.sort_values('Overall Ranking')
    
    # Prepare results with cleaned output
    results = {
        'target_solvent_properties': target.to_dict(),
        'by_similarity': top_similar[['Name','SMILES','Density', 'Dielectric', 'Dipole','Refractive Index', 'Melting point', 'Boiling point']],
        'by_environment': env_ranked[['Name', 'Environment Ranking']],
        'by_health': health_ranked[['Name', 'Health Ranking']],
        'by_safety': safety_ranked[['Name','Safety Ranking']],
        'by_overall_ranking': overall_ranked[['Name', 'Adjusted ranking','Overall Ranking']]
    }
    
    return results

result = rank_similar_solvents('O')
print(result['target_solvent_properties'])