"""Test configuration for pytest."""

import pytest
import pandas as pd
from unittest.mock import MagicMock, patch

@pytest.fixture
def mock_df():
    """Create a mock DataFrame for testing."""
    return pd.DataFrame([{
        'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1',
        'reactants': 'CC(=O)O.OC1=CC=CC=C1',
        'product': 'CC(=O)OC1=CC=CC=C1'
    }])

@pytest.fixture
def mock_solvent_df():
    """Create a mock solvent DataFrame for testing."""
    return pd.DataFrame({
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

# Add these new fixtures to replace retrogsf_test.py

@pytest.fixture
def mock_retrosynthesis():
    """Mock for retrosynthesis_reaction_smiles function."""
    with patch('retrogsf.retrogsf.AiZynthExpander') as mock_expander_class:
        # Configure the mock instance
        mock_expander = MagicMock()
        mock_expander_class.return_value = mock_expander
        
        # Configure the expansion_policy and filter_policy
        mock_expander.expansion_policy = MagicMock()
        mock_expander.filter_policy = MagicMock()
        
        # Configure do_expansion to return mock data
        mock_reactions = [[MagicMock()]]
        mock_reactions[0][0].metadata = {
            'mapped_reaction_smiles': 'CC(=O)OC1=CC=CC=C1>>CC(=O)O.OC1=CC=CC=C1',
            'reactants': 'CC(=O)O.OC1=CC=CC=C1',
            'product': 'CC(=O)OC1=CC=CC=C1'
        }
        mock_expander.do_expansion.return_value = mock_reactions
        
        yield

@pytest.fixture
def mock_reaction_info():
    """Mock for rxn_info function."""
    with patch('retrogsf.retrogsf.Reaction') as mock_reaction_class:
        # Configure the mock
        mock_reaction = MagicMock()
        mock_reaction_class.return_value = mock_reaction
        
        # Set up the get_reaction_info method to return a dictionary with NAME key
        mock_info = {'NAME': 'Ester Hydrolysis'}
        mock_reaction.get_reaction_info.return_value = mock_info
        
        yield


@pytest.fixture
def mock_similar_solvents():
    """Mock for rank_similar_solvents function."""
    with patch('retrogsf.retrogsf.pd.read_csv') as mock_read_csv:
        # Configure the mock with all required columns
        mock_df = pd.DataFrame({
            'Name': ['Water', 'Ethanol', 'Acetonitrile'],
            'SMILES': ['O', 'CCO', 'CC#N'],
            'Density': [1.0, 0.8, 0.7],
            'Dielectric': [80, 24, 36],
            'Dipole': [1.85, 1.7, 3.9],
            'Refractive Index': [1.33, 1.36, 1.34],
            'Melting point': [0, -114, -45],
            'Boiling point': [100, 78, 82],
            'Environment Ranking': [1, 2, 3],
            'Health Ranking': [1, 2, 3],
            'Safety Ranking': [1, 2, 3],
            'Overall Ranking': [1, 2, 3],
            'Adjusted ranking': ['Recommended', 'Recommended', 'Problematic']
        })
        mock_read_csv.return_value = mock_df
        
        yield
