import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from pathlib import Path

from retrogsf.retrogsf import (
    retrosynthesis_reaction_smiles,
    rxn_info,
    get_solvents_for_reaction,
    rank_similar_solvents
)


@pytest.fixture
def mock_expander_result():
    """Mock the result from AiZynthExpander"""
    metadata = [{
        'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1',
        'reactants': 'CC(=O)O.OC1=CC=CC=C1',
        'product': 'CC(=O)OC1=CC=CC=C1'
    }]
    return [[(MagicMock(metadata=item)) for item in metadata]]


@pytest.fixture
def mock_reaction_info():
    """Mock reaction info result"""
    return {"NAME": "Ester Hydrolysis", "CLASS": "Hydrolysis"}


@patch('retrogsf.retrogsf.AiZynthExpander')
def test_retrosynthesis_reaction_smiles(mock_expander_class, mock_expander_result):
    """Test retrosynthesis_reaction_smiles function"""
    # Setup mock
    mock_instance = mock_expander_class.return_value
    mock_instance.do_expansion.return_value = mock_expander_result
    
    # Call function
    result = retrosynthesis_reaction_smiles("CC(=O)OC1=CC=CC=C1")
    
    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert 'mapped_reaction_smiles' in result.columns
    assert result.iloc[0]['mapped_reaction_smiles'] == 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'


@patch('retrogsf.retrogsf.Reaction')
def test_rxn_info(mock_reaction_class, mock_reaction_info):
    """Test rxn_info function"""
    # Setup mock
    mock_instance = mock_reaction_class.return_value
    mock_instance.get_reaction_info.return_value = mock_reaction_info
    
    # Create test dataframe
    df = pd.DataFrame([{'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'}])
    
    # Call function
    result = rxn_info(df)
    
    # Assertions
    assert result == "Ester Hydrolysis"


@patch('retrogsf.retrogsf.genai.Client')
def test_get_solvents_for_reaction(mock_client_class):
    """Test get_solvents_for_reaction function with mocked Gemini API"""
    # Setup mock response
    mock_response = MagicMock()
    mock_response.text = "O, CCO, CC#N"
    
    # Setup mock client
    mock_client = mock_client_class.return_value
    mock_client.models.generate_content.return_value = mock_response
    
    # Call function
    result = get_solvents_for_reaction("Ester Hydrolysis")
    
    # Assertions
    assert result == "O, CCO, CC#N"
    mock_client.models.generate_content.assert_called_once()
    # Check that the model name is correct
    args, kwargs = mock_client.models.generate_content.call_args
    assert kwargs.get('model') == "gemini-2.0-flash"
    # Check that the prompt contains the reaction name
    assert "Ester Hydrolysis" in kwargs.get('contents', [''])[0]


@patch('retrogsf.retrogsf.pd.read_csv')
def test_rank_similar_solvents(mock_read_csv):
    """Test rank_similar_solvents function"""
    # Setup mock dataframe
    mock_df = pd.DataFrame({
        'SMILES': ['O', 'CCO', 'CC#N'],
        'Name': ['Water', 'Ethanol', 'Acetonitrile'],
        'Adjusted ranking': ['Recommended', 'Recommended', 'Problematic'],
        'Environment Ranking': [1, 2, 3],
        'Health Ranking': [1, 2, 3],
        'Safety Ranking': [1, 2, 3],
        'Density': [1.0, 0.8, 0.7],
        'Dielectric': [80, 24, 36],
        'Dipole': [1.8, 1.7, 3.9],
        'Refractive Index': [1.33, 1.36, 1.34],
        'Melting point': [0, -114, -45],
        'Boiling point': [100, 78, 82]
    })
    mock_read_csv.return_value = mock_df
    
    # Call function
    result = rank_similar_solvents('O')
    
    # Assertions
    assert isinstance(result, dict)
    assert 'by_similarity' in result
    assert 'by_environment' in result
    assert 'by_health' in result
    assert 'by_safety' in result
    assert 'by_overall_ranking' in result


def test_rank_similar_solvents_not_found():
    """Test rank_similar_solvents when SMILES not found"""
    # Call function with a non-existent SMILES
    result = rank_similar_solvents('NONEXISTENT')
    
    # Assertions
    assert isinstance(result, str)
    assert "Error: No solvent found with SMILES 'NONEXISTENT'" in result