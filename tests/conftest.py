"""Test configuration for pytest."""

import pytest
import pandas as pd
from unittest.mock import MagicMock

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