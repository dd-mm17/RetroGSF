import sys
from pathlib import Path

retrogsf_path = Path(__file__).resolve().parent.parent / 'src'
sys.path.append(str(retrogsf_path))

import streamlit as st
from rdkit.Chem.Draw import IPythonConsole
from rxn_insight.reaction import Reaction
from dotenv import load_dotenv
import os
import pandas as pd
import numpy as np
from aizynthfinder.aizynthfinder import AiZynthExpander
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import requests
from streamlit_ketcher import st_ketcher
from retrogsf import retrosynthesis_reaction_smiles, rxn_info, get_solvents_for_reaction, rank_similar_solvents, unmap_reaction_smiles

# ==== RXN Drawing ====
def draw_reaction_with_solvent(reactants, products, solvent_text, img_width=400):

    rxn_smiles = f"{reactants}>>{products}"
    rxn = Chem.rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    

    img = Draw.ReactionToImage(rxn, subImgSize=(img_width, int(img_width / 2)))


    img_with_text = Image.new("RGBA", (img.width, img.height + 40), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img_with_text)
    img_with_text.paste(img, (0, 40))


    try:
        font = ImageFont.truetype("arial.ttf", 16)
    except OSError:
        font = ImageFont.load_default()
        

    text_width, _ = draw.textsize(solvent_text, font=font)
    draw.text(((img.width - text_width) / 2, 10), solvent_text, fill="black", font=font)

    return img_with_text

# ==== APP STREAMLIT ====
st.title("🧪 RetroGSF")
config_input=st.text_input("Enter the path to your config.yml file:")
input_method = st.radio("Choose input method:", ["Enter SMILES", "Draw Molecule"])

if input_method == "Enter SMILES":
    smiles_input = st.text_input("Enter a product SMILES string:")
else:
    # Ketcher drawing interface
    with st.expander("Draw Molecule From SMILES (optional)", expanded=True):
        initial_smiles = st.text_input("**Initial SMILES for drawing**", "")
    ketcher_smiles = st_ketcher(initial_smiles, height=600)
    
    with st.expander("Generated SMILES from Drawing"):
        st.markdown(f"**SMILES:** {ketcher_smiles}")
    
    smiles_input = ketcher_smiles

data_path_input = st.text_input("Enter the path to the solvent data file (optional):")

if smiles_input:
    try:
        df = retrosynthesis_reaction_smiles(smiles_input, config_path=config_input)
        reaction_name = rxn_info(df)
        rxn_smiles=df.iloc[0]['mapped_reaction_smiles']
        unmapped_rxn = unmap_reaction_smiles(rxn_smiles)
        reactants, reagents, products = unmapped_rxn.split('>')

        # Seperate reactants
        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split('.')]
        #st.subheader(f"Reaction Name or Class: {reaction_name}")

        solvents = get_solvents_for_reaction(reaction_name)

        # get_iupac_name function
        def get_iupac_name(smiles):
            try:
                url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
                response = requests.get(url)
                if response.status_code == 200:
                    return response.text.strip()
                else:
                    return "Name not found"
            except Exception as e:
                return f"Error: {str(e)}"

        solvents_name = get_iupac_name(solvents).capitalize()

        st.subheader("🧬 Reaction with Suggested Solvent")
        img = draw_reaction_with_solvent(products, reactants, solvents) 
        st.image(img, caption=f"Suggested solvent: {solvents_name}", use_container_width=True)

        st.subheader("🌐 Informations about the retro-synthesis:")
        st.write(f"Reaction SMILES: {rxn_smiles}")
        st.write(f"Name or class of the reaction: {reaction_name}")

        # Pass the custom path to the function if provided
        results = rank_similar_solvents(solvents, data_path=data_path_input if data_path_input else None)
        if isinstance(results, str):
                st.error(results)
        else :
            tabs = st.tabs([
                 "🎯 Target solvent properties"
            ])     

            with tabs[0]:
                st.write("Top solvents by **target solvent properties**:")
                st.dataframe(results['target_solvent_properties'])

        st.subheader("🔬 Similar & Safer Solvent Recommendations")

        if isinstance(results, str):
            st.error(results)
        else :
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
