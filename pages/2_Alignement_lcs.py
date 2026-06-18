"""
Lot 2 : Alignement LCS (Plus Longue Sous-Séquence Commune)
"""
import streamlit as st
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.alignment import align_reads
from aff.ui import inject_styles, page_header, sidebar_brand, sidebar_section, lot_chart, LOT_COLOR

st.set_page_config(page_title="Lot 2 — Alignement LCS", layout="wide")
inject_styles(lot=2)
CHART = lot_chart(2)

page_header(
    "Lot 2 — Alignement LCS",
    "LCS (Plus Longue Sous-Séquence Commune) par programmation dynamique (DP).",
)

with st.sidebar:
    sidebar_brand(lot=2)
    sidebar_section("Rappels")
    st.markdown("""
    **LCS** : Plus Longue Sous-Séquence Commune.

    **DP** : Programmation Dynamique.

    **bp** : paire de bases.

    Symboles : `|` match, `.` mismatch, espace = gap.
    """)

input_mode = st.radio(
    "Saisie",
    ["Manuelle", "Fichier FASTA", "Exemple"],
    horizontal=True,
)

seq1_val, seq2_val = "", ""

if "Manuelle" in input_mode:
    col1, col2 = st.columns(2)
    with col1:
        seq1_val = st.text_area("Read 1", value="ATGGAGTTCAACACCAATGTCACAGGCGAATCG", height=80)
    with col2:
        seq2_val = st.text_area("Read 2", value="TCACAGGCGAATCGTACAACCAGCACATCATTG", height=80)

elif "FASTA" in input_mode:
    uploaded = st.file_uploader("FASTA (2 séquences min.)", type=["fasta", "fa", "txt"])
    if uploaded:
        seqs, current = [], ""
        for line in uploaded.read().decode('utf-8', errors='ignore').split('\n'):
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
            seq1_val, seq2_val = seqs[0][:200], seqs[1][:200]
        else:
            st.error("Au moins 2 séquences requises.")

else:
    example = st.selectbox("Exemple", [
        "Chevauchement court",
        "Chevauchement long",
        "Avec erreur",
        "Sans chevauchement",
    ])
    examples = {
        "Chevauchement court": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGT",
        ),
        "Chevauchement long": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACAT",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCA",
        ),
        "Avec erreur": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG".replace('T', 'C', 1),
        ),
        "Sans chevauchement": (
            "AAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
    }
    seq1_val, seq2_val = examples[example]

run_align = st.button(
    "Calculer l'alignement",
    type="primary",
    disabled=not (seq1_val.strip() and seq2_val.strip()),
)

if run_align or (seq1_val and seq2_val and "Exemple" in input_mode):
    s1, s2 = seq1_val.strip().upper(), seq2_val.strip().upper()
    if len(s1) < 5 or len(s2) < 5:
        st.error("Séquences de 5 bases minimum.")
    else:
        result = align_reads(s1, s2)
        viz = result['visual']

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Longueur read 1 (bp)", len(s1))
        c2.metric("Longueur read 2 (bp)", len(s2))
        c3.metric("LCS (bp)", result['lcs_length'])
        c4.metric("Identité LCS", f"{result['identity'] * 100:.1f} %")

        col_left, col_right = st.columns([3, 2])

        with col_left:
            st.markdown("#### Alignement")
            s1_al, s2_al, match_l = viz['seq1_aligned'], viz['seq2_aligned'], viz['match_line']
            for i in range(0, len(s1_al), 60):
                block = {'s1': s1_al[i:i + 60], 's2': s2_al[i:i + 60],
                         'match': match_l[i:i + 60], 'start': i}
                s1_html, s2_html = "", ""
                for c1c, c2c, mc in zip(block['s1'], block['s2'], block['match']):
                    col = LOT_COLOR[2] if mc == '|' else CHART["error"] if mc == '.' else CHART["muted"]
                    s1_html += f'<span style="color:{col}">{c1c}</span>'
                    s2_html += f'<span style="color:{col}">{c2c}</span>'
                match_html = ''.join(
                    f'<span style="color:{LOT_COLOR[2] if c == "|" else CHART["error"] if c == "." else CHART["muted"]}">{c}</span>'
                    for c in block['match']
                )
                pos = block['start']
                st.markdown(f"""
                <div class="align-block">
<span style="color:{LOT_COLOR[2]};font-weight:600">read1 {pos:>4}</span>  {s1_html}
<span style="color:{CHART['muted']}">           </span>  {match_html}
<span style="color:{LOT_COLOR[2]};font-weight:600">read2 {pos:>4}</span>  {s2_html}
                </div>
                """, unsafe_allow_html=True)

            st.caption(
                f"Chevauchement : {result['overlap_length']} bp "
                f"({result['overlap_identity'] * 100:.1f} % identité) · offset {viz['offset']}"
            )
            st.code(result['lcs'] or "(aucune)", language="text")

        with col_right:
            st.markdown("#### Matrice DP")
            dp = result['matrix']
            n_rows, n_cols = dp.shape
            MAX_DISP = 25
            dp_disp = dp[:min(n_rows, MAX_DISP), :min(n_cols, MAX_DISP)]

            fig, ax = plt.subplots(
                figsize=(min(n_cols, MAX_DISP) * 0.5, min(n_rows, MAX_DISP) * 0.45),
                facecolor=CHART['bg'],
            )
            ax.set_facecolor(CHART['face'])
            im = ax.imshow(dp_disp, cmap='Purples', aspect='auto')
            ax.set_xticks(range(dp_disp.shape[1]))
            ax.set_xticklabels(['-'] + list(s2[:MAX_DISP - 1]), color=CHART['text'], fontsize=7)
            ax.set_yticks(range(dp_disp.shape[0]))
            ax.set_yticklabels(['-'] + list(s1[:MAX_DISP - 1]), color=CHART['text'], fontsize=7)
            if min(n_rows, MAX_DISP) <= 15 and min(n_cols, MAX_DISP) <= 15:
                for i in range(dp_disp.shape[0]):
                    for j in range(dp_disp.shape[1]):
                        ax.text(j, i, str(dp_disp[i, j]), ha='center', va='center', fontsize=7,
                                color="#FFFFFF" if dp_disp[i, j] > dp_disp.max() * 0.55 else CHART["title"])
            ax.set_title(f"dp[{n_rows}×{n_cols}]", color=CHART['title'], fontsize=9)
            plt.colorbar(im, ax=ax).ax.tick_params(colors=CHART['text'])
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()
