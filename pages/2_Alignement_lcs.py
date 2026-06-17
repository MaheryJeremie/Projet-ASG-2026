"""
Lot 2 : Alignement par Programmation Dynamique (LCS)
"""
import streamlit as st
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.alignment import align_reads
from aff.ui import inject_styles, page_header, align_legend, CHART

st.set_page_config(page_title="Lot 2 — Alignement LCS", layout="wide")
inject_styles()

page_header(
    "Lot 2 — Alignement LCS",
    "Plus Longue Sous-Séquence Commune entre deux reads. Complexité O(n × m) en temps et en espace.",
)
st.divider()

with st.sidebar:
    st.markdown("### Rappel LCS")
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

st.markdown("### Saisie des séquences")

input_mode = st.radio(
    "Mode de saisie",
    ["Saisie manuelle", "Upload fichier FASTA", "Exemples prédéfinis"],
    horizontal=True,
)

seq1_val, seq2_val = "", ""

if "Saisie manuelle" in input_mode:
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Read 1**")
        seq1_val = st.text_area(
            "", value="ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            height=80, key="seq1_input", label_visibility="collapsed",
        )
    with col2:
        st.markdown("**Read 2**")
        seq2_val = st.text_area(
            "", value="TCACAGGCGAATCGTACAACCAGCACATCATTG",
            height=80, key="seq2_input", label_visibility="collapsed",
        )

elif "Upload" in input_mode:
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
            st.success(f"{len(seqs)} séquences chargées. Utilisation des 2 premières.")
        else:
            st.error("Le fichier doit contenir au moins 2 séquences.")

else:
    example = st.selectbox("Choisir un exemple", [
        "Chevauchement court (20 bp)",
        "Chevauchement long (50 bp)",
        "Séquences similaires avec erreurs",
        "Séquences non chevauchantes",
    ])

    examples = {
        "Chevauchement court (20 bp)": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGT",
        ),
        "Chevauchement long (50 bp)": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACAT",
            "CGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCA",
        ),
        "Séquences similaires avec erreurs": (
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG",
            "ATGGAGTTCAACACCAATGTCACAGGCGAATCG".replace('T', 'C', 1),
        ),
        "Séquences non chevauchantes": (
            "AAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
    }
    seq1_val, seq2_val = examples[example]

    col1, col2 = st.columns(2)
    col1.code(f"Read 1: {seq1_val}")
    col2.code(f"Read 2: {seq2_val}")

st.divider()

run_align = st.button(
    "Calculer l'alignement", type="primary",
    disabled=not (seq1_val.strip() and seq2_val.strip()),
)

if run_align or (seq1_val and seq2_val and "Exemples" in input_mode):
    s1 = seq1_val.strip().upper()
    s2 = seq2_val.strip().upper()

    if len(s1) < 5 or len(s2) < 5:
        st.error("Les séquences doivent comporter au moins 5 bases.")
    else:
        with st.spinner("Calcul de la matrice LCS..."):
            result = align_reads(s1, s2)

        st.markdown("### Résultats d'alignement")

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Read 1", f"{len(s1)} bp")
        c2.metric("Read 2", f"{len(s2)} bp")
        c3.metric("LCS", f"{result['lcs_length']} bp")
        c4.metric("Identité LCS", f"{result['identity'] * 100:.1f} %")
        c5.metric(
            "Chevauchement", f"{result['overlap_length']} bp",
            delta=f"{result['overlap_identity'] * 100:.1f} % id.",
        )

        st.divider()

        col_left, col_right = st.columns([3, 2])

        with col_left:
            st.markdown("#### Alignement visuel (chevauchement)")

            viz = result['visual']
            s1_al = viz['seq1_aligned']
            s2_al = viz['seq2_aligned']
            match_l = viz['match_line']

            block_size = 60
            blocks = []
            for i in range(0, len(s1_al), block_size):
                blocks.append({
                    's1': s1_al[i:i + block_size],
                    's2': s2_al[i:i + block_size],
                    'match': match_l[i:i + block_size],
                    'start': i,
                })

            for b in blocks:
                s1_html = ""
                s2_html = ""
                for c1c, c2c, mc in zip(b['s1'], b['s2'], b['match']):
                    if mc == '|':
                        s1_html += f'<span style="color:{CHART["success"]}">{c1c}</span>'
                        s2_html += f'<span style="color:{CHART["success"]}">{c2c}</span>'
                    elif mc == '.':
                        s1_html += f'<span style="color:{CHART["error"]}">{c1c}</span>'
                        s2_html += f'<span style="color:{CHART["error"]}">{c2c}</span>'
                    else:
                        s1_html += f'<span style="color:{CHART["muted"]}">{c1c}</span>'
                        s2_html += f'<span style="color:{CHART["muted"]}">{c2c}</span>'

                pos = b['start']
                match_html = ''.join(
                    f'<span style="color:{CHART["success"] if c == "|" else CHART["error"] if c == "." else CHART["muted"]}">{c}</span>'
                    for c in b['match']
                )
                st.markdown(f"""
                <div class="align-block">
<span style="color:{CHART['accent']}">Read1 {pos:>4}</span>  {s1_html}
<span style="color:{CHART['muted']}">          </span>  {match_html}
<span style="color:{CHART['accent']}">Read2 {pos:>4}</span>  {s2_html}
                </div>
                """, unsafe_allow_html=True)

            st.caption(align_legend(viz['offset'], viz['identity']))

            st.markdown("#### Plus Longue Sous-Séquence Commune (LCS)")
            st.code(result['lcs'] if result['lcs'] else "(aucune)", language="text")

        with col_right:
            st.markdown("#### Matrice LCS (DP)")

            dp = result['matrix']
            n_rows, n_cols = dp.shape

            MAX_DISP = 25
            dp_disp = dp[:min(n_rows, MAX_DISP), :min(n_cols, MAX_DISP)]

            fig, ax = plt.subplots(
                figsize=(min(n_cols, MAX_DISP) * 0.5, min(n_rows, MAX_DISP) * 0.45),
                facecolor=CHART['bg'],
            )
            ax.set_facecolor(CHART['bg'])

            im = ax.imshow(dp_disp, cmap='YlOrRd', aspect='auto')

            ax.set_xticks(range(dp_disp.shape[1]))
            ax.set_xticklabels(['-'] + list(s2[:MAX_DISP - 1]), color=CHART['text'], fontsize=7)
            ax.set_yticks(range(dp_disp.shape[0]))
            ax.set_yticklabels(['-'] + list(s1[:MAX_DISP - 1]), color=CHART['text'], fontsize=7)

            if min(n_rows, MAX_DISP) <= 15 and min(n_cols, MAX_DISP) <= 15:
                for i in range(dp_disp.shape[0]):
                    for j in range(dp_disp.shape[1]):
                        ax.text(
                            j, i, str(dp_disp[i, j]),
                            ha='center', va='center', fontsize=7,
                            color=CHART['bg'] if dp_disp[i, j] > dp_disp.max() * 0.5 else CHART['legend_text'],
                        )

            ax.set_title(
                f"dp[{n_rows}×{n_cols}]" + (" (tronquée)" if n_rows > MAX_DISP or n_cols > MAX_DISP else ""),
                color=CHART['title'], fontsize=9,
            )
            plt.colorbar(im, ax=ax).ax.tick_params(colors=CHART['text'])
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()

            st.caption(
                f"Matrice complète : {result['dp_dimensions']}  \n"
                f"Complexité spatiale : O({len(s1)} × {len(s2)}) = O({len(s1) * len(s2):,}) cellules"
            )

            st.markdown("#### Récapitulatif")
            summary = {
                "Longueur Read 1": f"{len(s1)} bp",
                "Longueur Read 2": f"{len(s2)} bp",
                "Score LCS": result['lcs_length'],
                "Identité globale": f"{result['identity'] * 100:.1f} %",
                "Longueur chevauchement": f"{result['overlap_length']} bp",
                "Identité chevauchement": f"{result['overlap_identity'] * 100:.1f} %",
                "Offset optimal": viz['offset'],
            }
            for label, val in summary.items():
                col_l, col_v = st.columns([2, 1])
                col_l.caption(label)
                col_v.markdown(f"**{val}**")
