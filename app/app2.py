import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import io
from retrogsf import retrosynthesis_reaction_smiles, rxn_info, get_solvents_for_reaction, rank_similar_solvents

# ==== TES FONCTIONS A IMPORTER ====
# from ton_module import retrosynthesis_reaction_smiles, get_solvents_for_reaction, rank_similar_solvents

# ==== STUBS POUR TESTS ====
def retrosynthesis_reaction_smiles(smiles: str, config_path: str = "config.yml") -> pd.DataFrame:
    return pd.DataFrame({
        "reactants": ["CCO.CCBr"],   # exemple : éthanol + bromoéthane
        "products": [smiles],
        "reaction_name": ["SN2_substitution"]
    })

def get_solvents_for_reaction(rxn_name):
    return ["Ethanol", "Water"]

def rank_similar_solvents(target_smiles, data_path='RetroGSF/data/Solvant_properties_with_smile.csv', n_recommendations=5):
    return [("Methanol", 0.91), ("Acetone", 0.89), ("DMSO", 0.88)]

# ==== FONCTION POUR AFFICHER LA RÉACTION ====
def draw_reaction_with_solvent(reactants, products, solvent_text):
    rxn_smiles = f"{reactants}>>{products}"
    rxn = Chem.rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=(300, 150))
    
    # Ajouter le nom du solvant au-dessus de la flèche
    img_with_text = Image.new("RGBA", (img.width, img.height + 40), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img_with_text)
    img_with_text.paste(img, (0, 40))

    # Centrer le solvant
    font = ImageFont.load_default()
    text_width, _ = draw.textsize(solvent_text, font=font)
    draw.text(((img.width - text_width) / 2, 10), solvent_text, fill="black", font=font)

    return img_with_text

# ==== INTERFACE STREAMLIT ====
st.title("🧪 Retrosynthesis + Solvent Suggestion")

smiles_input = st.text_input("🔬 Enter a product SMILES string:")

if smiles_input:
    try:
        # Étape 1 : Retrosynthèse
        df = retrosynthesis_reaction_smiles(smiles_input)
        reactants = df.iloc[0]["reactants"]
        products = df.iloc[0]["products"]
        rxn_name = df.iloc[0]["reaction_name"]

        # Étape 2 : Solvant suggéré
        solvents = get_solvents_for_reaction(rxn_name)
        solvent = solvents[0] if solvents else "Unknown"

        # Étape 3 : Dessiner la réaction
        reaction_image = draw_reaction_with_solvent(reactants, products, solvent)
        st.image(reaction_image, caption=f"Reaction using {solvent}")

        # Étape 4 : Tableau de solvants similaires
        ranked_solvents = rank_similar_solvents(smiles_input)
        df_ranked = pd.DataFrame(ranked_solvents, columns=["Solvent", "Similarity"])
        st.subheader("🔍 Similar Solvent Recommendations")
        st.dataframe(df_ranked)

    except Exception as e:
        st.error(f"❌ Error: {e}")

