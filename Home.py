"""
ASG-2026 : Pipeline d'Assemblage De Novo
"""
import streamlit as st

from aff.ui import inject_styles, page_header, module_card, sidebar_brand, sidebar_section

st.set_page_config(
    page_title="ASG-2026 | Pipeline d'Assemblage De Novo",
    layout="wide",
    initial_sidebar_state="expanded",
)

inject_styles()

with st.sidebar:
    sidebar_brand()
    sidebar_section("Navigation")
    st.caption("Choisir un lot dans le menu.")

page_header(
    "ASG-2026 — Assemblage de novo",
    "Trois lots : ingestion, alignement LCS, assemblage par filtre de Bloom.",
)

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown(
        module_card(
            1, "Ingestion et qualité",
            "Lecture FASTQ, histogramme des k-mers, conversion FASTA.",
            ["FASTQ", "k-mers", "Phred"],
            lot=1,
        ),
        unsafe_allow_html=True,
    )

with col2:
    st.markdown(
        module_card(
            2, "Alignement LCS",
            "Plus Longue Sous-Séquence Commune, matrice de programmation dynamique.",
            ["LCS", "DP", "chevauchement"],
            lot=2,
        ),
        unsafe_allow_html=True,
    )

with col3:
    st.markdown(
        module_card(
            3, "Assemblage Bloom",
            "Filtre de Bloom, graphe de de Bruijn implicite, contigs.",
            ["Bloom", "de Bruijn", "contigs"],
            lot=3,
        ),
        unsafe_allow_html=True,
    )

st.code("""
FASTQ  →  Lot 1 (k-mers, qualité)  →  Lot 2 (LCS)  →  Lot 3 (Bloom, contigs)
""", language="text")
