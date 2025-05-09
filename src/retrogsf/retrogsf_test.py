"""
Test-friendly version of retrogsf module with minimal dependencies.
This module is used only for testing purposes.
"""

import pandas as pd

def mock_retrosynthesis_reaction_smiles(smiles: str, config_path: str = "config.yml") -> pd.DataFrame:
    """Mock version of retrosynthesis_reaction_smiles for testing"""
    return pd.DataFrame([{
        'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1',
        'reactants': 'CC(=O)O.OC1=CC=CC=C1',
        'product': 'CC(=O)OC1=CC=CC=C1'
    }])

def mock_rxn_info(df: pd.DataFrame) -> str:
    """Mock version of rxn_info for testing"""
    return "Ester Hydrolysis"

def mock_get_solvents_for_reaction(rxn_name: str) -> str:
    """Mock version of get_solvents_for_reaction for testing"""
    return "O, CCO, CC#N"

def mock_rank_similar_solvents(target_smiles: str, data_path: str = 'SHE_data_with_smiles.csv', n_recommendations: int = 5):
    """Mock version of rank_similar_solvents for testing"""
    if target_smiles == 'NONEXISTENT':
        return f"Error: No solvent found with SMILES '{target_smiles}' in the dataset."
        
    # Return a mock result
    return {
        'target_solvent_properties': {'Name': 'Water', 'SMILES': 'O'},
        'by_similarity': pd.DataFrame({
            'Name': ['Ethanol', 'Acetonitrile'],
            'SMILES': ['CCO', 'CC#N'],
            'Density': [0.8, 0.7],
            'Dielectric': [24, 36],
            'Dipole': [1.7, 3.9],
            'Refractive Index': [1.36, 1.34],
            'Melting point': [-114, -45],
            'Boiling point': [78, 82]
        }),
        'by_environment': pd.DataFrame({
            'Name': ['Ethanol', 'Acetonitrile'],
            'Environment Ranking': [2, 3]
        }),
        'by_health': pd.DataFrame({
            'Name': ['Ethanol', 'Acetonitrile'],
            'Health Ranking': [2, 3]
        }),
        'by_safety': pd.DataFrame({
            'Name': ['Ethanol', 'Acetonitrile'],
            'Safety Ranking': [2, 3]
        }),
        'by_overall_ranking': pd.DataFrame({
            'Name': ['Ethanol', 'Acetonitrile'],
            'Adjusted ranking': ['Recommended', 'Problematic'],
            'Overall Ranking': [2, 3]
        })
    }