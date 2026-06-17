"""
ASG-2026 : Pipeline d'Assemblage De Novo
Application principale Streamlit
"""
import streamlit as st

from aff.ui import inject_styles, page_header, module_card

st.set_page_config(
    page_title="ASG-2026 | Pipeline d'Assemblage De Novo",
    layout="wide",
    initial_sidebar_state="expanded",
)

inject_styles()

page_header(
    "ASG-2026 — Pipeline d'assemblage de novo",
    "Reconstruction de séquences génomiques à partir de reads courts. "
    "Trois étapes : ingestion et contrôle qualité, alignement LCS, assemblage par graphe de de Bruijn implicite.",
)
st.divider()

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown(
        module_card(
            1,
            "Ingestion et qualité",
            "Lecture FASTQ/FASTA, comptage des k-mers, histogramme de fréquence, estimation du taux d'erreur.",
            ["FASTQ → FASTA", "k-mers", "Phred"],
        ),
        unsafe_allow_html=True,
    )

with col2:
    st.markdown(
        module_card(
            2,
            "Alignement LCS",
            "Plus Longue Sous-Séquence Commune entre deux reads, visualisation du chevauchement.",
            ["LCS O(n²)", "chevauchement", "matrice DP"],
        ),
        unsafe_allow_html=True,
    )

with col3:
    st.markdown(
        module_card(
            3,
            "Assemblage Bloom",
            "Filtre de Bloom et traversée du graphe de de Bruijn implicite (approche Minia 2).",
            ["Bloom filter", "de Bruijn", "contigs"],
        ),
        unsafe_allow_html=True,
    )

st.divider()

st.markdown("#### Architecture du pipeline")
st.code("""
Fichiers FASTQ/FASTA
        │
        ▼
┌───────────────────┐
│  Lot 1 : Ingestion │  → k-mer counting → histogramme fréquence → seuil "solide"
│  & Contrôle Qualité│
└─────────┬─────────┘
          │  séquences filtrées
          ▼
┌───────────────────┐
│  Lot 2 : Alignement│  → matrice LCS O(n×m) → score + position chevauchement
│  Programmation Dyn.│
└─────────┬─────────┘
          │  validation locale
          ▼
┌───────────────────┐
│  Lot 3 : Assemblage│  → Bloom Filter (k-mers solides)
│  De Bruijn Implicite│  → Traversée on-the-fly → CONTIGS
└───────────────────┘
""", language="text")

st.info("Les trois modules sont accessibles dans le menu latéral.")

with st.expander("Jeu de test (toy dataset)"):
    st.markdown("""
    Chaque module inclut un jeu de test synthétique : séquence de référence d'environ 500 bp,
    10 000 reads simulés, taux d'erreur 1 %.

    **Validation** : la séquence reconstruite doit atteindre une identité ≥ 98 % par rapport à la cible.

    Référence : fragment synthétique du gène Spike d'un coronavirus fictif.
    """)
