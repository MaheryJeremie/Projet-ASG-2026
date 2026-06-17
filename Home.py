"""
ASG-2026 : Pipeline d'Assemblage De Novo
Application principale Streamlit
"""
import streamlit as st

st.set_page_config(
    page_title="ASG-2026 | Pipeline d'Assemblage De Novo",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── CSS personnalisé ───────────────────────────────────────────────────────
st.markdown("""
<style>
    /* Palette : fond sombre lab + accents bio-vert */
    [data-testid="stAppViewContainer"] {
        background: #0d1117;
        color: #e6edf3;
    }
    [data-testid="stSidebar"] {
        background: #161b22;
        border-right: 1px solid #30363d;
    }
    h1, h2, h3 { color: #58a6ff; }
    .lot-card {
        background: #161b22;
        border: 1px solid #30363d;
        border-left: 4px solid #3fb950;
        border-radius: 8px;
        padding: 1.2rem 1.5rem;
        margin-bottom: 1rem;
    }
    .lot-card h3 { color: #3fb950; margin: 0 0 0.4rem 0; font-size: 1rem; }
    .lot-card p { color: #8b949e; margin: 0; font-size: 0.9rem; }
    .badge {
        display: inline-block;
        background: #21262d;
        border: 1px solid #3fb950;
        color: #3fb950;
        font-size: 0.75rem;
        padding: 2px 8px;
        border-radius: 12px;
        margin-right: 6px;
    }
    .metric-row {
        display: flex; gap: 1rem; margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ─── Header ─────────────────────────────────────────────────────────────────
st.markdown("## 🧬 ASG-2026")
st.markdown("### Pipeline d'Assemblage De Novo")
st.markdown(
    "Suite logicielle de reconstruction de séquences génomiques à partir de reads courts. "
    "Implémente un pipeline complet : ingestion → alignement → assemblage par graphe de de Bruijn implicite."
)
st.divider()

# ─── Description des lots ───────────────────────────────────────────────────
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    <div class="lot-card">
        <h3>📦 Lot 1 — Ingestion & Qualité</h3>
        <p>Lecture de fichiers FASTQ/FASTA, analyse des k-mers et histogramme de fréquence pour estimer le taux d'erreur.</p>
        <br/>
        <span class="badge">FASTQ → FASTA</span>
        <span class="badge">k-mers</span>
        <span class="badge">Qualité Phred</span>
    </div>
    """, unsafe_allow_html=True)

with col2:
    st.markdown("""
    <div class="lot-card">
        <h3>🔗 Lot 2 — Alignement LCS</h3>
        <p>Programmation dynamique pour calculer la Plus Longue Sous-Séquence Commune entre deux reads et visualiser le chevauchement.</p>
        <br/>
        <span class="badge">LCS O(n²)</span>
        <span class="badge">Chevauchement</span>
        <span class="badge">Matrice DP</span>
    </div>
    """, unsafe_allow_html=True)

with col3:
    st.markdown("""
    <div class="lot-card">
        <h3>⚙️ Lot 3 — Assemblage Bloom</h3>
        <p>Moteur d'assemblage memory-efficient via Filtre de Bloom et traversée on-the-fly du graphe de de Bruijn implicite.</p>
        <br/>
        <span class="badge">Bloom Filter</span>
        <span class="badge">de Bruijn</span>
        <span class="badge">Minia 2</span>
    </div>
    """, unsafe_allow_html=True)

st.divider()

# ─── Architecture technique ─────────────────────────────────────────────────
st.markdown("#### 🏗️ Architecture du Pipeline")
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

st.info("👈 **Naviguez** entre les modules via le menu latéral.", icon="ℹ️")

# ─── Toy Dataset info ───────────────────────────────────────────────────────
with st.expander("🧪 Toy Dataset de validation"):
    st.markdown("""
    Un jeu de test synthétique est intégré dans chaque module. Il utilise une séquence
    de référence courte (~500 bp) avec **10 000 reads simulés** et un taux d'erreur de **1%**.
    
    **Critère de validation** : La séquence reconstruite doit avoir une identité ≥ 98% 
    avec la séquence cible.
    
    La séquence de référence utilisée est un fragment synthétique du gène de la 
    protéine de spicule (Spike) d'un coronavirus fictif.
    """)