"""
Lot 2 : Alignement par Programmation Dynamique (LCS)
- Saisie ou upload de deux reads
- Calcul de la LCS
- Affichage lisible de l'alignement
"""
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.alignment import align_reads, lcs_matrix

st.set_page_config(page_title="Lot 2 — Alignement LCS", page_icon="🔗", layout="wide")

st.markdown("""
<style>
[data-testid="stAppViewContainer"] { background: #0d1117; color: #e6edf3; }
[data-testid="stSidebar"] { background: #161b22; border-right: 1px solid #30363d; }
h1, h2, h3 { color: #58a6ff; }
.align-block {
    font-family: 'Courier New', monospace;
    font-size: 0.85rem;
    background: #161b22;
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 1rem 1.2rem;
    line-height: 1.8;
    overflow-x: auto;
    white-space: pre;
}
.match { color: #3fb950; }
.mismatch { color: #f85149; }
.gap { color: #8b949e; }
</style>
""", unsafe_allow_html=True)

st.markdown("## 🔗 Lot 2 — Alignement par Programmation Dynamique")
st.markdown(
    "Calcule la **Plus Longue Sous-Séquence Commune (LCS)** entre deux reads "
    "et visualise leur chevauchement. Complexité : **O(n × m)** en temps et espace."
)
st.divider()

# ─── Sidebar : théorie ───────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### 📖 Rappel LCS")
    st.markdown(r"""
    **Récurrence :**
    
    Si $s_1[i] = s_2[j]$ :
    $$dp[i][j] = dp[i-1][j-1] + 1$$
    
    Sinon :
    $$dp[i][j] = \max(dp[i-1][j],\ dp[i][j-1])$$
    
    **Score d'identité :**
    $$\text{Id} = \frac{|LCS|}{\max(|s_1|, |s_2|)}$$
    
    **Complexité :**
    - Temps : $O(n \cdot m)$
    - Espace : $O(n \cdot m)$
    """)
    
    st.divider()
    st.markdown("### Symboles d'alignement")
    st.markdown("""
    `|` → Match parfait  
    `.` → Mismatch  
    ` ` → Gap  
    """)

# ─── Saisie des reads ────────────────────────────────────────────────────────
st.markdown("### 📥 Saisie des séquences")

input_mode = st.radio("Mode de saisie", ["✏️ Saisie manuelle", "📁 Upload fichier FASTA",
                                          "🧪 Exemples prédéfinis"], horizontal=True)

seq1_val, seq2_val = "", ""

if "✏️ Saisie" in input_mode:
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Read 1**")
        seq1_val = st.text_area("", value="ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
                                 height=80, key="seq1_input", label_visibility="collapsed")
    with col2:
        st.markdown("**Read 2**")
        seq2_val = st.text_area("", value="TCACAGGCGAATCGTACAACCAGCACATCATTG",
                                 height=80, key="seq2_input", label_visibility="collapsed")

elif "📁 Upload" in input_mode:
    uploaded = st.file_uploader("Fichier FASTA (au moins 2 séquences)", type=["fasta", "fa", "txt"])
    if uploaded:
        content = uploaded.read().decode('utf-8', errors='ignore')
        seqs = []
        current = ""
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current:
                    seqs.append(current)
                current = ""
            elif line:
                current += line.upper()
        if current:
            seqs.append(current)
        
        if len(seqs) >= 2:
            seq1_val = seqs[0][:200]
            seq2_val = seqs[1][:200]
            st.success(f"✅ {len(seqs)} séquences chargées. Utilisation des 2 premières.")
        else:
            st.error("Le fichier doit contenir au moins 2 séquences.")

else:  # Exemples prédéfinis
    example = st.selectbox("Choisir un exemple", [
        "Chevauchement court (20 bp)",
        "Chevauchement long (50 bp)",
        "Séquences similaires avec erreurs",
        "Séquences non chevauchantes",
    ])
    
    examples = {
        "Chevauchement court (20 bp)": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGT"
        ),
        "Chevauchement long (50 bp)": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACAT",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCA"
        ),
        "Séquences similaires avec erreurs": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG".replace('T','C', 1)
        ),
        "Séquences non chevauchantes": (
            "AAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTTTTTTTTTTT"
        ),
    }
    seq1_val, seq2_val = examples[example]
    
    col1, col2 = st.columns(2)
    col1.code(f"Read 1: {seq1_val}")
    col2.code(f"Read 2: {seq2_val}")

# ─── Bouton d'alignement ────────────────────────────────────────────────────
st.divider()

run_align = st.button("⚡ Calculer l'alignement", type="primary",
                       disabled=not (seq1_val.strip() and seq2_val.strip()))

if run_align or (seq1_val and seq2_val and "🧪" in input_mode):
    s1 = seq1_val.strip().upper()
    s2 = seq2_val.strip().upper()
    
    if len(s1) < 5 or len(s2) < 5:
        st.error("Les séquences doivent comporter au moins 5 bases.")
    else:
        with st.spinner("Calcul de la matrice LCS..."):
            result = align_reads(s1, s2)
        
        # ─── Métriques ───────────────────────────────────────────────────────
        st.markdown("### 📊 Résultats d'alignement")
        
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Read 1", f"{len(s1)} bp")
        c2.metric("Read 2", f"{len(s2)} bp")
        c3.metric("LCS", f"{result['lcs_length']} bp")
        c4.metric("Identité LCS", f"{result['identity']*100:.1f}%")
        c5.metric("Chevauchement", f"{result['overlap_length']} bp",
                  delta=f"{result['overlap_identity']*100:.1f}% id.")
        
        st.divider()
        
        col_left, col_right = st.columns([3, 2])
        
        with col_left:
            # ─── Alignement visuel ──────────────────────────────────────────
            st.markdown("#### 🔍 Alignement visuel (chevauchement)")
            
            viz = result['visual']
            s1_al = viz['seq1_aligned']
            s2_al = viz['seq2_aligned']
            match_l = viz['match_line']
            
            # Découper en blocs de 60
            block_size = 60
            blocks = []
            for i in range(0, len(s1_al), block_size):
                blocks.append({
                    's1': s1_al[i:i+block_size],
                    's2': s2_al[i:i+block_size],
                    'match': match_l[i:i+block_size],
                    'start': i,
                })
            
            for b in blocks:
                # Colorisation HTML
                s1_html = ""
                s2_html = ""
                for c1c, c2c, mc in zip(b['s1'], b['s2'], b['match']):
                    if mc == '|':
                        s1_html += f'<span style="color:#3fb950">{c1c}</span>'
                        s2_html += f'<span style="color:#3fb950">{c2c}</span>'
                    elif mc == '.':
                        s1_html += f'<span style="color:#f85149">{c1c}</span>'
                        s2_html += f'<span style="color:#f85149">{c2c}</span>'
                    else:
                        s1_html += f'<span style="color:#8b949e">{c1c}</span>'
                        s2_html += f'<span style="color:#8b949e">{c2c}</span>'
                
                pos = b['start']
                st.markdown(f"""
                <div class="align-block">
<span style="color:#58a6ff">Read1 {pos:>4}</span>  {s1_html}
<span style="color:#444c56">          </span>  {''.join('<span style="color:'+('#3fb950' if c=='|' else '#f85149' if c=='.' else '#444c56')+'">'+c+'</span>' for c in b['match'])}
<span style="color:#58a6ff">Read2 {pos:>4}</span>  {s2_html}
                </div>
                """, unsafe_allow_html=True)
            
            # Légende
            st.markdown(
                "🟢 Match &nbsp;|&nbsp; 🔴 Mismatch &nbsp;|&nbsp; "
                f"⬛ Gap &nbsp;|&nbsp; **Offset** : {viz['offset']} "
                f"&nbsp;|&nbsp; **Identité de chevauchement** : {viz['identity']*100:.1f}%"
            )
            
            # LCS
            st.markdown("#### 🧬 Plus Longue Sous-Séquence Commune (LCS)")
            st.code(result['lcs'] if result['lcs'] else "(aucune)", language="text")
        
        with col_right:
            # ─── Matrice DP ─────────────────────────────────────────────────
            st.markdown("#### 🗂️ Matrice LCS (DP)")
            
            dp = result['matrix']
            n_rows, n_cols = dp.shape
            
            # Limiter l'affichage à 25×25 pour lisibilité
            MAX_DISP = 25
            dp_disp = dp[:min(n_rows, MAX_DISP), :min(n_cols, MAX_DISP)]
            
            fig, ax = plt.subplots(
                figsize=(min(n_cols, MAX_DISP) * 0.5, min(n_rows, MAX_DISP) * 0.45),
                facecolor='#0d1117'
            )
            ax.set_facecolor('#0d1117')
            
            im = ax.imshow(dp_disp, cmap='YlOrRd', aspect='auto')
            
            # Labels des axes
            ax.set_xticks(range(dp_disp.shape[1]))
            ax.set_xticklabels(['-'] + list(s2[:MAX_DISP-1]),
                                color='#8b949e', fontsize=7)
            ax.set_yticks(range(dp_disp.shape[0]))
            ax.set_yticklabels(['-'] + list(s1[:MAX_DISP-1]),
                                color='#8b949e', fontsize=7)
            
            # Valeurs dans les cellules
            if min(n_rows, MAX_DISP) <= 15 and min(n_cols, MAX_DISP) <= 15:
                for i in range(dp_disp.shape[0]):
                    for j in range(dp_disp.shape[1]):
                        ax.text(j, i, str(dp_disp[i, j]),
                                ha='center', va='center', fontsize=7,
                                color='#0d1117' if dp_disp[i, j] > dp_disp.max() * 0.5 else '#e6edf3')
            
            ax.set_title(
                f"dp[{n_rows}×{n_cols}]" + (" (tronquée)" if n_rows > MAX_DISP or n_cols > MAX_DISP else ""),
                color='#58a6ff', fontsize=9
            )
            plt.colorbar(im, ax=ax).ax.tick_params(colors='#8b949e')
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()
            
            st.caption(
                f"Matrice complète : {result['dp_dimensions']}  \n"
                f"Complexité spatiale : O({len(s1)} × {len(s2)}) = O({len(s1)*len(s2):,}) cellules"
            )
            
            # Récapitulatif
            st.markdown("#### 📋 Récapitulatif")
            summary = {
                "Longueur Read 1": f"{len(s1)} bp",
                "Longueur Read 2": f"{len(s2)} bp",
                "Score LCS": result['lcs_length'],
                "Identité globale": f"{result['identity']*100:.1f}%",
                "Longueur chevauchement": f"{result['overlap_length']} bp",
                "Identité chevauchement": f"{result['overlap_identity']*100:.1f}%",
                "Offset optimal": viz['offset'],
            }
            for label, val in summary.items():
                col_l, col_v = st.columns([2, 1])
                col_l.caption(label)
                col_v.markdown(f"**{val}**")