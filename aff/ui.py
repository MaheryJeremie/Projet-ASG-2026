"""Couche présentation Streamlit — styles et helpers d'affichage."""
import streamlit as st

CHART = {
    "bg": "#1a1a1c",
    "face": "#242428",
    "accent": "#6b9ec2",
    "success": "#7dba8a",
    "error": "#c97a7a",
    "warning": "#c9a85c",
    "text": "#8a8a94",
    "title": "#d4d4d4",
    "grid": "#2e2e34",
    "spine": "#3a3a40",
    "legend_bg": "#242428",
    "legend_edge": "#3a3a40",
    "legend_text": "#d4d4d4",
    "muted": "#5c5c66",
}

_CSS = """
<style>
:root {
    --bg: #1a1a1c;
    --surface: #242428;
    --border: #3a3a40;
    --text: #d4d4d4;
    --muted: #8a8a94;
    --accent: #6b9ec2;
    --mono: Consolas, 'Courier New', monospace;
}

[data-testid="stAppViewContainer"] {
    background: var(--bg);
    color: var(--text);
}

[data-testid="stSidebar"] {
    background: var(--surface);
    border-right: 1px solid var(--border);
}

.block-container {
    max-width: 1100px;
    padding-top: 1.5rem;
}

h1, h2, h3 { color: var(--text) !important; font-weight: 600 !important; }
h4, h5, h6 { color: var(--muted) !important; }
p, li { color: var(--muted); line-height: 1.6; }

[data-testid="stSidebar"] p,
[data-testid="stSidebar"] li {
    font-size: 0.9rem;
}

.katex-display {
    padding: 0.5rem 0;
    overflow-x: auto;
}

.module-card {
    border: 1px solid var(--border);
    border-left: 3px solid var(--accent);
    padding: 1rem 1.1rem;
    margin-bottom: 0.25rem;
    height: 100%;
}
.module-card-title {
    color: var(--text);
    font-weight: 600;
    font-size: 0.95rem;
    margin-bottom: 0.45rem;
}
.module-n {
    color: var(--accent);
    font-variant-numeric: tabular-nums;
    margin-right: 0.35rem;
}
.module-card p {
    color: var(--muted);
    font-size: 0.875rem;
    margin: 0 0 0.65rem 0;
    line-height: 1.5;
}
.module-tags {
    color: var(--muted);
    font-family: var(--mono);
    font-size: 0.75rem;
}

.align-block {
    font-family: var(--mono);
    font-size: 0.8rem;
    background: var(--surface);
    border: 1px solid var(--border);
    padding: 0.85rem 1rem;
    line-height: 1.85;
    overflow-x: auto;
    white-space: pre;
    margin-bottom: 0.4rem;
}

.contig-block {
    font-family: var(--mono);
    font-size: 0.78rem;
    background: var(--surface);
    border: 1px solid var(--border);
    padding: 0.7rem 0.9rem;
    margin-bottom: 0.4rem;
    line-height: 1.5;
    word-break: break-all;
}

.stat-line {
    padding: 0.35rem 0;
    border-bottom: 1px solid var(--border);
    font-size: 0.875rem;
    color: var(--muted);
}
.stat-line:last-child { border-bottom: none; }
.stat-line b { color: var(--text); font-weight: 500; }

[data-testid="stCode"] pre {
    font-family: var(--mono) !important;
    font-size: 0.8rem !important;
}

hr { border-color: var(--border) !important; }
</style>
"""


def inject_styles() -> None:
    st.markdown(_CSS, unsafe_allow_html=True)


def page_header(title: str, subtitle: str = "") -> None:
    st.markdown(f"## {title}")
    if subtitle:
        st.caption(subtitle)


def module_card(number: int, title: str, description: str, tags: list[str]) -> str:
    tags_html = " · ".join(tags)
    return f"""
    <div class="module-card">
        <div class="module-card-title"><span class="module-n">{number}.</span>{title}</div>
        <p>{description}</p>
        <div class="module-tags">{tags_html}</div>
    </div>
    """


def stat_line(label: str, value: str) -> str:
    return f'<div class="stat-line">{label} : <b>{value}</b></div>'


def align_legend(offset: int, identity: float) -> str:
    return (
        f"Match en vert · Mismatch en rouge · Gap en gris — "
        f"Offset {offset} · Identité de chevauchement {identity * 100:.1f} %"
    )
