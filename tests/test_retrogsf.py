import pytest
import pandas as pd
from unittest.mock import patch

# Use the test-friendly versions for testing
from retrogsf.retrogsf_test import (
    mock_retrosynthesis_reaction_smiles,
    mock_rxn_info,
    mock_get_solvents_for_reaction,
    mock_rank_similar_solvents
)

# Simple test to verify the package can be imported
def test_import():
    """Test that the package can be imported"""
    import retrogsf
    assert retrogsf is not None

# Test the mock functions directly
def test_retrosynthesis_reaction_smiles():
    """Test mock_retrosynthesis_reaction_smiles function"""
    result = mock_retrosynthesis_reaction_smiles("CC(=O)OC1=CC=CC=C1")
    
    assert isinstance(result, pd.DataFrame)
    assert 'mapped_reaction_smiles' in result.columns
    assert result.iloc[0]['mapped_reaction_smiles'] == 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'

def test_rxn_info():
    """Test mock_rxn_info function"""
    df = pd.DataFrame([{'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'}])
    result = mock_rxn_info(df)
    
    assert result == "Ester Hydrolysis"

def test_get_solvents_for_reaction():
    """Test mock_get_solvents_for_reaction function"""
    result = mock_get_solvents_for_reaction("Ester Hydrolysis")
    
    assert result == "O, CCO, CC#N"

def test_rank_similar_solvents():
    """Test mock_rank_similar_solvents function"""
    result = mock_rank_similar_solvents('O')
    
    assert isinstance(result, dict)
    assert 'by_similarity' in result
    assert 'by_environment' in result
    assert 'by_health' in result
    assert 'by_safety' in result
    assert 'by_overall_ranking' in result

def test_rank_similar_solvents_not_found():
    """Test mock_rank_similar_solvents when SMILES not found"""
    result = mock_rank_similar_solvents('NONEXISTENT')
    
    assert isinstance(result, str)
    assert "Error: No solvent found with SMILES 'NONEXISTENT'" in result
