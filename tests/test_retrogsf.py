import pytest
import pandas as pd
import re

# Import the actual functions to test
from retrogsf.retrogsf import (
    retrosynthesis_reaction_smiles,
    rxn_info,
    get_solvents_for_reaction,
    rank_similar_solvents,
    unmap_reaction_smiles
)

# Test unmap_reaction_smiles (no external dependencies, can test directly)
def test_unmap_reaction_smiles():
    """Test unmap_reaction_smiles function"""
    mapped = "C[C:1](=[O:2])[O:3][C:4]1=[C:5][C:6]=[C:7][C:8]=[C:9]1>>[C:1][C:2](=[O:3])[O:4].[O:5][C:6]1=[C:7][C:8]=[C:9][C:10]=[C:11]1"
    result = unmap_reaction_smiles(mapped)
    
    # Check that all mapping numbers are removed
    assert not re.search(r':\d+\]', result)
    # Check that the basic structure is preserved
    assert ">>" in result
    assert result.count(">") == 2

# Test retrosynthesis_reaction_smiles with mocked AiZynthExpander
def test_retrosynthesis_reaction_smiles(mock_retrosynthesis):
    """Test retrosynthesis_reaction_smiles function with mocked AiZynthExpander"""
    test_smiles = "CC(=O)OC1=CC=CC=C1"
    
    # Call the function
    result = retrosynthesis_reaction_smiles(test_smiles, config_path="dummy_path")
    
    # Verify the result
    assert isinstance(result, pd.DataFrame)
    assert 'mapped_reaction_smiles' in result.columns
    assert result.iloc[0]['mapped_reaction_smiles'] == 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'

# Test rxn_info with a mock for Reaction class
def test_rxn_info(mock_reaction_info, mock_df):
    """Test rxn_info function with mocked Reaction"""
    # Call the function
    result = rxn_info(mock_df)
    
    # Verify the result
    assert result == "Ester Hydrolysis"

# Test get_solvents_for_reaction with mocked Google Gemini API
def test_get_solvents_for_reaction(mock_solvent_recommendation):
    """Test get_solvents_for_reaction function with mocked Gemini API"""
    # Call the function
    result = get_solvents_for_reaction("Ester Hydrolysis")
    
    # Verify the result
    assert result == "O, CCO, CC#N"

def test_rank_similar_solvents(mock_similar_solvents):
    """Test rank_similar_solvents function with mocked pandas operations"""
    # Call the function with a valid SMILES
    result = rank_similar_solvents('O')
    
    # Verify the result is a dictionary with expected keys
    assert isinstance(result, dict)
    assert 'target_solvent_properties' in result
    assert 'by_similarity' in result
    assert 'by_environment' in result
    assert 'by_health' in result
    assert 'by_safety' in result
    assert 'by_overall_ranking' in result
