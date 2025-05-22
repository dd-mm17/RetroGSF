import pytest
import pandas as pd
import re
from unittest.mock import patch, MagicMock

from retrogsf.retrogsf import (
    retrosynthesis_reaction_smiles,
    rxn_info,
    get_solvents_for_reaction,
    rank_similar_solvents,
    unmap_reaction_smiles
)

def test_unmap_reaction_smiles():
    """Test unmap_reaction_smiles function"""
    mapped = "C[C:1](=[O:2])[O:3][C:4]1=[C:5][C:6]=[C:7][C:8]=[C:9]1>>[C:1][C:2](=[O:3])[O:4].[O:5][C:6]1=[C:7][C:8]=[C:9][C:10]=[C:11]1"
    result = unmap_reaction_smiles(mapped)
    assert not re.search(r':\d+\]', result)
    assert ">>" in result
    assert result.count(">") == 2

# Test retrosynthesis_reaction_smiles with mocked AiZynthExpander
def test_retrosynthesis_reaction_smiles(mock_retrosynthesis):
    """Test retrosynthesis_reaction_smiles function with mocked AiZynthExpander"""
    test_smiles = "CC(=O)OC1=CC=CC=C1"
    result = retrosynthesis_reaction_smiles(test_smiles, config_path="dummy_path")
    assert isinstance(result, pd.DataFrame)
    assert 'mapped_reaction_smiles' in result.columns
    assert result.iloc[0]['mapped_reaction_smiles'] == 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1'

# Test rxn_info with a mock for Reaction class
def test_rxn_info(mock_reaction_info, mock_df):
    """Test rxn_info function with mocked Reaction"""
    result = rxn_info(mock_df)
    assert result == "Ester Hydrolysis"

# Test get_solvents_for_reaction with mocked Google Gemini API
def test_get_solvents_for_reaction():
    """Test get_solvents_for_reaction function with mocked Gemini API"""
    with patch("google.generativeai.GenerativeModel") as MockModel:
        mock_model_instance = MockModel.return_value
        mock_response = MagicMock()
        mock_response.text = "O"
        mock_model_instance.generate_content.return_value = mock_response

        # Patch genai.configure to do nothing
        with patch("google.generativeai.configure"):
            # Patch os.environ.get to return a fake API key
            with patch("os.environ.get", return_value="fake_api_key"):
                result = get_solvents_for_reaction("Ester Hydrolysis")
                assert result == "O"

def test_rank_similar_solvents(mock_similar_solvents):
    """Test rank_similar_solvents function with mocked pandas operations"""
    result = rank_similar_solvents('O')
    if isinstance(result, str):
        assert result == "Warning: No safe solvents found within range"
    else:
        assert isinstance(result, dict)
        assert 'target_solvent_properties' in result
        assert 'by_similarity' in result
        assert 'by_environment' in result
        assert 'by_health' in result
        assert 'by_safety' in result
        assert 'by_overall_ranking' in result
