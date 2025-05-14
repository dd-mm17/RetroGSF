import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
from retrogsf import retrosynthesis_reaction_smiles, rxn_info, get_solvents_for_reaction, rank_similar_solvents

def retrosynthesis_reaction_smiles(smiles: str, config_path: str = "config.yml") -> pd.DataFrame:
    return pd.DataFrame({
        "reactants": ["CCO.CCBr"],
        "products": [smiles],
        "reaction_name": ["SN2_substitution"]
    })

def get_solvents_for_reaction(rxn_name):
    return ["Ethanol", "Water"]

def rank_similar_solvents(target_smiles, data_path='SHE_data_with_smiles.csv', n_recommendations=5):
    data = {
        'Name': ['Methanol', 'Acetone', 'DMSO'],
        'SMILES': ['CO', 'CC(=O)C', 'CS(=O)C'],
        'Density': [0.79, 0.79, 1.1],
        'Dielectric': [33, 21, 47],
        'Dipole': [1.7, 2.88, 4.0],
        'Refractive Index': [1.328, 1.358, 1.479],
        'Melting point': [-97, -95, 18],
        'Boiling point': [65, 56, 189],
        'Environment Ranking': [2, 3, 4],
        'Health Ranking': [2, 3, 3],
        'Safety Ranking': [1, 2, 3],
        'Adjusted ranking': ['Recommended', 'Recommended', 'Recommended'],
    }
    df = pd.DataFrame(data)
    return {
        'target_solvent_properties': {'Name': 'Ethanol', 'Density': 0.789, 'Boiling point': 78},
        'by_similarity': df[['Name','SMILES','Density', 'Dielectric', 'Dipole','Refractive Index', 'Melting point', 'Boiling point']],
        'by_environment': df[['Name', 'Environment Ranking']],
        'by_health': df[['Name', 'Health Ranking']],
        'by_safety': df[['Name','Safety Ranking']],
        'by_overall_ranking': df[['Name', 'Adjusted ranking','Environment Ranking']]
    }

# ==== RXN Drawing ====
def draw_reaction_with_solvent(reactants, products, solvent_text):
    rxn_smiles = f"{reactants}>>{products}"
    rxn = Chem.rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=(300, 150))

    img_with_text = Image.new("RGBA", (img.width, img.height + 40), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img_with_text)
    img_with_text.paste(img, (0, 40))

    font = ImageFont.load_default()
    text_width, _ = draw.textsize(solvent_text, font=font)
    draw.text(((img.width - text_width) / 2, 10), solvent_text, fill="black", font=font)

    return img_with_text

# ==== APP STREAMLIT ====
st.title("🧪 Retrosynthesis and Safer Solvent Recommendations")

smiles_input = st.text_input("Enter a product SMILES string:")

if smiles_input:
    try:
        # 1. Retrosynthèse
        df = retrosynthesis_reaction_smiles(smiles_input)
        reactants = df.iloc[0]["reactants"]
        products = df.iloc[0]["products"]
        rxn_name = df.iloc[0]["reaction_name"]

        # 2. Solvant principal proposé
        solvents = get_solvents_for_reaction(rxn_name)
        main_solvent = solvents[0] if solvents else "Unknown"

        # 3. Dessin de la réaction avec solvant
        st.subheader("🧪 Reaction with Suggested Solvent")
        img = draw_reaction_with_solvent(reactants, products, main_solvent)
        st.image(img, caption=f"Suggested solvent: {main_solvent}")

        # 4. Appel à rank_similar_solvents
        st.subheader("🔬 Similar & Safer Solvent Recommendations")
        results = rank_similar_solvents(main_solvent)

        if isinstance(results, str):
            st.error(results)
        else:
            # Affichage sous forme d'onglets pour chaque DataFrame
            tabs = st.tabs([
                "🔍 Similarity-based",
                "🌱 Environment Ranking",
                "🩺 Health Ranking",
                "⚠️ Safety Ranking",
                "⭐ Overall Ranking"
            ])

            with tabs[0]:
                st.write("Top solvents by **physical similarity**:")
                st.dataframe(results['by_similarity'])

            with tabs[1]:
                st.write("Top solvents by **environmental impact**:")
                st.dataframe(results['by_environment'])

            with tabs[2]:
                st.write("Top solvents by **health impact**:")
                st.dataframe(results['by_health'])

            with tabs[3]:
                st.write("Top solvents by **safety score**:")
                st.dataframe(results['by_safety'])

            with tabs[4]:
                st.write("Top solvents by **overall ranking**:")
                st.dataframe(results['by_overall_ranking'])

    except Exception as e:
        st.error(f"❌ Error: {e}")