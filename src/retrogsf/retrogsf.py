import os
import numpy as np
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rxn_insight.reaction import Reaction
from aizynthfinder.aizynthfinder import AiZynthExpander
from dotenv import load_dotenv
from google import genai
from google.genai import types
load_dotenv()

def retrosynthesis_reaction_smiles(smiles: str, config_path: str = "config.yml") -> pd.DataFrame:
    """
    Perform retrosynthesis and return a table of forward-ordered one-step Reaction SMILES.

    Args:
        smiles (str): Target molecule in SMILES format.
        config_path (str): Path to AiZynthFinder's config.yml file. (see their github for more info)

    Returns:
        pd.DataFrame: Table with step number, reactants, product, and Reaction SMILES.
    """

    p = Path("/Users/diego/Desktop/EPFL/Prog. in Chem/data_download/config.yml") # Change path to config file on git
    expander = AiZynthExpander(configfile=p)
    expander.expansion_policy.select("uspto")
    expander.filter_policy.select("uspto")
    reactions = expander.do_expansion(smiles)
    metadata = []
    for reaction_tuple in reactions:
        for reaction in reaction_tuple:
            metadata.append(reaction.metadata)
    df = pd.DataFrame(metadata)
    return df

def rxn_info (df: pd.DataFrame) -> str:
    """
    Get the reaction name or class from the reaction SMILES. 
    If name is 'OtherReaction', returns the class.

    Args:
        df (pd.DataFrame): DataFrame containing the reaction SMILES.

    Returns:
        str: Reaction name or class.
    """
    rxn_smiles=df.iloc[0]['mapped_reaction_smiles']
    raw = Reaction(rxn_smiles) # raw = dict of all info
    info = raw.get_reaction_info()
    if info.get("NAME") != "OtherReaction":
        name_class = info.get("NAME", "Unknown")
    else:
        name_class = info.get("CLASS", "Unknown")
    return name_class

def get_solvents_for_reaction(rxn_name):
    """
    Get recommended solvents for a given reaction name/type.
    
    Args:
        rxn_name (str): The name or type of the reaction
        
    Returns:
        str: Comma-separated SMILES strings of recommended solvents
    """
    known_solvents = [
        'O',
        'C1CC2(C(=O)CC1O2)O',
        'C(C(CO)O)O',
        'OCCO',
        'CC(O)CO' ,
        'OCCOCCO',
        'OCCN(CCO)CCO',
        'CS(=O)C',
        'CCO',
        'CC(C)O',
        'CCCC(=O)O',
        'CC#N',
        'OCCOCCOCCO',
        'C[N+](=O)[O-]',
        'CC(C)(C)O',
        'CCC(C)O',
        'CC1COC(=O)O1',
        'CN(C)P(=O)(N(C)C)N(C)C',
        'CC(C)CCO',
        'CC(=O)C',
        'COCCOCO',
        'COC(C)CO',
        'CCCCCO',
        'COC(C)CO',
        'CCCOCCO',
        'CC(=O)CC(C)(C)O',
        'OC1CCCCC1',
        'OCc1ccccc1',
        'CCC(C)(C)O',
        'c1ccncc1',
        'CCCC(C)O',
        'CCOCCOCCO',
        'CCC(CC)O',
        'CCC(C)CO',
        'OCCOCCOCCOCCO',
        'CC(=O)OC',
        'O=C1CCCC1',
        'CCC(=O)C',
        'CC(C)CC(C)O',
        'CCCCOCCO',
        'CN1CCOCC1',
        'COC(OC)OC',
        'CCCC(=O)C',
        'CCCCOCCOCCO',
        'CCC(=O)OC',
        'COC(C)C(=O)OC',
        'CCCOC=O',
        'CC(=O)OCC',
        'COC(=O)OC',
        'O=C1CCCCC1',
        'CCCCCCCO',
        'CCC(=O)CC',
        'CCCOC(C)CO',
        'CC(=O)OCCOC(=O)C',
        'CC(=O)OC(C)C',
        'CCCOC(C)C(=O)OC',
        'CCCOC(=O)C',
        'CCCCOC=O',
        'CC(=O)OCC(OC(=O)C)COC(=O)C',
        'CC(=O)CC(C)C',
        'CCOCC',
        'CCOCCOCCOC',
        'CC(Cl)Cl',
        'COC(C)(C)C',
        'C1CCSC1',
        'CCCCOCCO',
        'CC(=O)c1ccccc1',
        'CC(=O)OCC(C)C',
        'CC1CCOC1',
        'COCC(C)OCC(C)OCC(C)O',
        'CCCCOC(=O)C',
        'CCOC(OCC)OCC',
        'c1ccsc1',
        'CC(C)CCOC(=O)C',
        'CCCCOCCO(C=O)C',
        'CCCCCOC(=O)C',
        'CC(C)OC(C)C',
        'COC(C)(C)CC',
        'COc1ccccc1',
        'c1ccc(Cl)cc1',
        'CCCCCCOC(=O)C',
        'CCCCCOC(=O)CC',
        'c1ccc(Br)cc1',
        'CCCCOCCOCCO(C=O)C',
        'c1ccc(Cl)c(Cl)c1',
        'CC1CCCC1',
        'CCCOC(=O)CC',
        'CCc1ccccc1',
        'CC(C)CCC',
        'CCCCC',
        'C1CCCCC1',
        'Cc1ccc(C)cc1',
        'CC1CCCCC1',
        'CC(C)c1ccccc1',
        'CCC1CCCCC1',
        'FC1=C(F)C(F)=C(F)C(F)=C1C(F)(F)F',
        'Cc1cc(C)cc(C)c1',
        'CCCCc1ccccc1',
        'CC1=CCC(CC1)C(=C)C',
        'CCCCCCC',
        'CCCCCCCC',
        'C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2',
        'c1ccc(Oc2ccccc2)cc1',
        'CCCCCCCCC',
        'CCCCCCCCCC',
        'C1(F)(F)C2(F)(F)C(F)(F)C(F)(F)C1(F)C1(F)C(F)(F)C(F)(F)C2(F)C1(F)F',
        'C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C',
        'C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C',
        'C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1',
        'CCCO',
        'O=S1CCCC1',
        'c1ccnnc1',
        'CCCCO',
        'CC(C)c1ccc(C)cc1',
        'CC(C)CO',
        'O=C1CCCCN1C',
        'c1cc[nH]c1',
        'CC1(C)OCC(CO)O1',
        'c1cncnc1',
        'CN1CCCC(=O)N(C)C1',
        'CN1C(=O)C=CC1=O',
        'CC[N+](=O)[O-]',
        'CN1CCCCCC1=O',
        'CC1CCC(=O)O1',
        'CC(O)C(=O)OCC',
        'CCCCCCOH',
        'CC(=O)CC(=O)C',
        'CCOC=O',
        'CCCC#N',
        'Nc1ccccc1',
        'CCOc1ccccc1',
        'O=S1OCCO1',
        'C1CCN=C2CCCCN2CC1',
        'COC1COC2C1OCC2OC',
        'COC(C)OC',
        'C1COCCO1',
        'COS(=O)OC',
        'CC(C)C(=O)C',
        'CC1CCC(CO)CC1',
        'CC(=O)CC(=O)OCC',
        'CC(OC)COCC(C)O',
        'CCC[N+](=O)[O-]',
        'CC1(C)OCCO1',
        'Cc1ccccn1',
        'CCCCOCCOCCOCCO',
        'C1CCOCC1',
        'CCCCC(=O)C',
        'CC1=CC(=O)CC(C)(C)C1',
        'CC(=O)C(C)(C)C',
        'c1ccccc1C#N',
        'CC(OC)COC(=O)C',
        'ClC=CCl',
        'CC1CCCCC1=O',
        'CCC(=O)OCC',
        'Cc1cccc(C)n1',
        'CCCC(=O)OC',
        'CCCCCC#N',
        'COC(C)OC',
        'CC(OC(=O)C)COC(=O)C',
        'COC(=O)CCCCC(=O)OC',
        'CCCCC(CC)CO',
        'ClCC(Cl)Cl',
        'CCCOC(C)COC(C)CO',
        'COC(OC)(OC)OC',
        'CC(COC)OC',
        'COC(C)COC(C)COC',
        'CCCCCCOCCOCCOH',
        'CCCC(=O)OCC',
        'CC(C)C(=O)OCC(C)(C)CC(C)O',
        'CCOC(=O)OCC',
        'CNc1ccccc1',
        'CC(C)C(=O)C(C)C',
        'C1CSCS1',
        'CCOCCCC(=O)OCC',
        'ClC=CCl',
        'CCCCCCCCOH',
        'CC1CCCO1',
        'Fc1ccccc1',
        'CCCCOCCOCCO',
        'CCCOCCC',
        'COC(=O)c1ccccc1',
        'CC(OC)COC(C)C(=O)OC',
        'CCNc1ccccc1',
        'ClC(Cl)C(Cl)(Cl)Cl',
        'FC(F)(F)c1ccccc1',
        'C1CCCC1',
        'CN(C)c1ccccc1',
        'Ic1ccccc1',
        'CCCCC(=O)CCCC',
        'CC(C)CC(=O)CC(C)C',
        'Cc1ccccc1',
        'CCOc1ccccc1',
        'ClC(=C(Cl)Cl)Cl',
        'CCCCOCCCC',
        'Clc1cccc(Cl)c1',
        'CCCCOCCOCCOCCO',
        'CCN(C(C)C)C(C)C',
        'c1ccc2CCCc2c1',
        'FC(F)(F)Oc1ccccc1',
        'FC1=C(F)C(F)=C(F)C(F)=C1F',
        'CCN(CC)S(=O)(=O)N(CC)CC',
        'c1ccc2c(c1)CCCC2',
        'CC1=CC=C(C=C1)C(C)C',
        'C[Si](C)(C)O[Si](C)(C)C',
        'C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1',
        'NC=O',
        'FC(F)(F)C(=O)O',
        'O=CO',
        'CNC=O',
        'NCCN',
        'CO',
        'NCCO',
        'CC(=O)NC',
        'CC(=O)O',
        'OC(C(F)(F)F)C(F)(F)F',
        'CCC(=O)O',
        'OCC(F)(F)F',
        '__loader__CN(C)C(=O)N(C)C',
        'CC(C)C(=O)O',
        'COCCO',
        'CN(C)C=O',
        'CC(=O)N(C)C',
        'CC(C)N',
        'O=S1(=O)CCCC1',
        'C1COCCN1',
        'O=C1CCCN(C)1',
        'CCCN',
        'CCOCCO',
        'CCCC#N',
        'C1CCNC1',
        'C1OCOCO1',
        'COC=O',
        'CC(C)CN',
        'CC(C)(C)N',
        'CCCCN',
        'O=CC1=CC=CO1',
        'N#CCCCCCC#N',
        'CCCCCN',
        'CC(=O)OC(=O)C',
        'C1CCCO1',
        'CN(C)P(=O)(N(C)C)N(C)C',
        'C1CCNCC1',
        'COCCOC',
        'NCc1ccccc1',
        'NC1CCCCC1',
        'COCCOCCOC',
        'COCCOCCOCCOC',
        'CCNCC',
        'ClCCl',
        'CC(C)[N+](=O)[O-]',
        'CCOCCOC(=O)C',
        'ClC(Cl)Cl',
        'ClCCCl',
        'O=[N+]([O-])c1ccccc1',
        'c1ccc2ncccc2c1',
        'ClC(Cl)C(Cl)Cl',
        'ICI',
        'ClC(=C)Cl',
        'c1ccoc1',
        'CCCCNCCC',
        'CC(C)NC(C)C',
        'S=C=S',
        'CCCCNCCCC',
        'CCN(CC)CC',
        'CC(Cl)(Cl)Cl',
        'ClC(Cl)=CCl',
        'c1ccccc1',
        'O=C(OC(=O)C(F)(F)F)C(F)(F)F',
        'ClC(Cl)(Cl)Cl',
        'CCCCCC',
        'C1CCC2CCCCC2C1',
        'CCCCN(CCCC)CCCC']

    client = genai.Client(api_key=os.environ.get("GEMINI_API_KEY")) 
    
    prompt=f""" 
    1. Main goal and context: 
    You are an expert in assigning solvents to reactions and know which solvent can be used for a given reactant. 
    You are part of a retro-synthesis program which will be used in mostly in solvent prediction (this is your job). 
    After your answer, the rest of the code will evaluate the "greeness" (how sustainable the possible solvents are).
    I hence need you to output only the Simplified Molecular Input Line Entry System (SMILES) of up to three possible 
    solvents which can be used for the given reaction name (or class if the name cannot be extracted for you, 
    of course, the class is quite general therefore you can be more general for your answer too). 

    The solvent you must find is to do with the following reaction name/type: {rxn_name}

    2. Constraints and examples
    The solvent you propose must be part of this list: {known_solvents}

    If you cannot find two solvents, one will do. 
    You must output at least one solvent and everything you output must be in the list 
    and in smiles format!
    You MUST ONLY output the smiles of the solvents in the following format: 
    "solventsmiles_1, solventsmiles_2, solventsmiles_3"

    This is an example output for you to visualise with the SMILES: "CN(C)C=O, ClCCl, CS(C)=O"

    3. Problems
    If you are given a reaction name or type which you do now know how to answer, you MUST simply reply with "nan"
    """

    response = client.models.generate_content(model="gemini-2.0-flash", contents=[prompt])
    response_stripped = response.text.strip()
    return response_stripped

def rank_similar_solvents(target_smiles, data_path='SHE_data_with_smiles.csv', n_recommendations=5):
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
    target = df[df['Name'] == target_smiles].iloc[0]
    
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
    
    # Filter out hazardous and highly hazardous compounds
    df_filtered = df_working[~df_working['Adjusted ranking'].isin(['Hazardous', 'Highly Hazardous'])].copy()
    
    # Strict filtering for all ranking scores (must be ≤ 5)
    mask = (
        (df_filtered['Environment Ranking'] <= 5) &
        (df_filtered['Health Ranking'] <= 5) &
        (df_filtered['Safety Ranking'] <= 5)
    )
    df_filtered = df_filtered[mask].copy()
    
    # Check if we still have solvents after filtering
    if len(df_filtered) == 0:
        return "Error: No solvents found meeting safety and environmental criteria (scores ≤ 5)."
    
    # Check the physical state and create temperature range filters
    target_melting = target['Melting point']
    target_boiling = target['Boiling point']
    
    # Important: Use df_filtered instead of df for temperature compatibility check
    df_filtered['temp_range_compatible'] = (
        (df_filtered['Melting point'] <= target_melting + 25) &  # Melting point not more than 25°C higher
        (df_filtered['Boiling point'] >= target_boiling - 25)    # Boiling point not more than 25°C lower
    )
    
    # Filter based on temperature criteria
    df_filtered = df_filtered[df_filtered['temp_range_compatible']].copy()
    
    # If no solvents match the temperature criteria, inform the user
    if len(df_filtered) == 0:
        return f"Warning: No safe solvents found within ±25°C range of '{target_solvent}' (MP: {target_melting}°C, BP: {target_boiling}°C)"
    
    # Calculate similarity scores based on physical properties
    # We'll compute a weighted Euclidean distance for the properties we care about
    
    # Define the properties to compare and their weights
    # It is also possible to add more properties here like 

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
    
    # Convert ranking columns to numeric using .loc
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