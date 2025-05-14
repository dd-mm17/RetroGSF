import streamlit as st
import pandas as pd
from retrogsf import retrosynthesis_reaction_smiles, rank_similar_solvents

# Titre de l'application
st.title("🧪 RetroGSF - Outil de rétrosynthèse et substitution de solvants")

# Saisie de l'utilisateur
smiles_input = st.text_input("🔍 Entrez le SMILES du produit cible", "")

if smiles_input:
    st.subheader("🧬 Résultat de la rétrosynthèse")
    try:
        # Étape 1 : rétrosynthèse
        df_retro = retrosynthesis_reaction_smiles(smiles_input)
        st.write("Voici les réactions en une étape générées :")
        st.dataframe(df_retro[['reactants', 'product', 'mapped_reaction_smiles']])
        
        # Étape 2 : sélection d'un solvant cible dans les réactifs
        st.subheader("🧴 Recherche de solvants similaires")
        unique_reactants = df_retro['reactants'].unique().tolist()
        selected_smiles = st.selectbox("Choisissez un SMILES de solvant cible parmi les réactifs :", unique_reactants)
        
        if selected_smiles:
            # Étape 3 : recommandation de solvants
            results = rank_similar_solvents(selected_smiles)
            if isinstance(results, str):
                st.error(results)
            else:
                st.write("✔️ Solvant cible :")
                st.json(results['target_solvent_properties'])

                st.write("📊 Solvants similaires classés par similarité :")
                st.dataframe(results['by_similarity'])

                st.write("🌱 Classement par critères :")
                col1, col2 = st.columns(2)
                with col1:
                    st.write("🟢 Environnement")
                    st.dataframe(results['by_environment'])
                with col2:
                    st.write("❤️ Santé")
                    st.dataframe(results['by_health'])

                st.write("🔥 Sécurité")
                st.dataframe(results['by_safety'])

                st.write("⭐ Classement global")
                st.dataframe(results['by_overall_ranking'])
    except Exception as e:
        st.error(f"Une erreur est survenue : {e}")