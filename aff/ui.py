"""
Présentation Streamlit — thème clair unifié, fort contraste, webdesign 2026.
"""
from __future__ import annotations

import streamlit as st

# Accents saturés (contraste élevé sur fond blanc)
LOT_COLOR = {1: "#2563EB", 2: "#7C3AED", 3: "#059669"}
LOT_SOFT = {1: "#DBEAFE", 2: "#EDE9FE", 3: "#D1FAE5"}
LOT_LABEL = {
    1: "Ingestion & qualité",
    2: "Alignement LCS",
    3: "Assemblage Bloom",
}

CHART = {
    "bg": "#FAFAFA",
    "face": "#FFFFFF",
    "accent": "#2563EB",
    "success": "#059669",
    "error": "#DC2626",
    "warning": "#D97706",
    "text": "#3F3F46",
    "title": "#09090B",
    "grid": "#E4E4E7",
    "spine": "#A1A1AA",
    "legend_bg": "#FFFFFF",
    "legend_edge": "#D4D4D8",
    "legend_text": "#18181B",
    "muted": "#71717A",
}

_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap');

:root {
    --bg: #FAFAFA;
    --surface: #FFFFFF;
    --border: #D4D4D8;
    --border-strong: #A1A1AA;
    --text: #09090B;
    --muted: #3F3F46;
    --dim: #71717A;
    --accent: #2563EB;
    --accent-soft: #DBEAFE;
    --sans: 'Inter', system-ui, sans-serif;
    --mono: 'JetBrains Mono', Consolas, monospace;
    --radius: 10px;
    --shadow: 0 1px 2px rgba(9,9,11,0.05);
}

html, body, [class*="css"] { font-family: var(--sans); }

/* ── Shell unifié (main + sidebar clairs) ─────────────────────────────── */
[data-testid="stAppViewContainer"] {
    background: var(--bg);
    color: var(--text);
}

[data-testid="stHeader"] { background: var(--surface); }

.block-container {
    max-width: 1140px;
    padding-top: 1.75rem;
    padding-bottom: 3rem;
}

h1, h2, h3 {
    color: var(--text) !important;
    font-weight: 700 !important;
    letter-spacing: -0.025em;
}
h4, h5, h6 { color: var(--muted) !important; font-weight: 600 !important; }
p, li { color: var(--muted); line-height: 1.65; font-size: 0.925rem; }

.katex-display {
    background: var(--surface);
    border: 2px solid var(--border);
    border-left: 3px solid var(--accent);
    border-radius: var(--radius);
    padding: 0.75rem 1rem;
    margin: 0.5rem 0;
    overflow-x: auto;
    color: var(--text);
}

/* ── En-têtes ──────────────────────────────────────────────────────────── */
.page-header {
    margin-bottom: 1.5rem;
    padding-bottom: 1rem;
    border-bottom: 2px solid var(--border);
}
.page-header h2 {
    font-size: 1.55rem !important;
    margin: 0 0 0.3rem !important;
    color: var(--text) !important;
}
.page-header .sub {
    color: var(--muted);
    font-size: 0.925rem;
    font-weight: 500;
    margin: 0;
}

/* ── Cartes ──────────────────────────────────────────────────────────── */
.dash-card {
    background: var(--surface);
    border: 2px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem 1.35rem;
    height: 100%;
    box-shadow: var(--shadow);
}
.dash-card:hover {
    border-color: var(--card-accent, var(--accent));
}
.dash-card-top { margin-bottom: 0.7rem; }
.dash-card-num {
    font-family: var(--mono);
    font-size: 0.68rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    color: var(--card-accent, var(--accent));
    background: var(--card-soft, var(--accent-soft));
    padding: 0.22rem 0.55rem;
    border-radius: 5px;
    border: 1px solid color-mix(in srgb, var(--card-accent, var(--accent)) 25%, transparent);
}
.dash-card h3 {
    font-size: 1.02rem !important;
    margin: 0 0 0.4rem !important;
    color: var(--text) !important;
}
.dash-card p {
    font-size: 0.875rem;
    font-weight: 500;
    margin: 0 0 0.75rem;
    line-height: 1.5;
    color: var(--muted);
}
.dash-card-tags {
    font-family: var(--mono);
    font-size: 0.7rem;
    font-weight: 600;
    color: var(--dim);
}

.panel {
    background: var(--surface);
    border: 2px solid var(--border);
    border-radius: var(--radius);
    padding: 1.25rem 1.35rem;
    margin-bottom: 1rem;
    box-shadow: var(--shadow);
}
.panel-title {
    font-size: 0.95rem;
    font-weight: 700;
    color: var(--text);
    margin: 0 0 0.75rem;
}

.align-block, .contig-block {
    font-family: var(--mono);
    font-size: 0.78rem;
    font-weight: 500;
    background: #F4F4F5;
    border: 2px solid var(--border);
    border-radius: 8px;
    padding: 0.85rem 1rem;
    line-height: 1.85;
    overflow-x: auto;
    white-space: pre;
    margin-bottom: 0.5rem;
    color: var(--text);
}
.contig-block {
    white-space: normal;
    word-break: break-all;
    line-height: 1.5;
}

.stat-line {
    display: flex;
    justify-content: space-between;
    padding: 0.5rem 0;
    border-bottom: 1px solid var(--border);
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--muted);
}
.stat-line:last-child { border-bottom: none; }
.stat-line b {
    font-family: var(--mono);
    font-size: 0.82rem;
    color: var(--text);
    font-weight: 700;
}

/* ── Sidebar claire ──────────────────────────────────────────────────── */
[data-testid="stSidebar"] {
    background: var(--surface) !important;
    border-right: 2px solid var(--border) !important;
}

[data-testid="stSidebar"] > div:first-child {
    background: var(--surface);
}

[data-testid="stSidebar"] .block-container {
    padding-top: 1.5rem;
}

[data-testid="stSidebar"] p,
[data-testid="stSidebar"] li {
    color: var(--muted) !important;
    font-weight: 500;
    font-size: 0.875rem;
}

[data-testid="stSidebar"] .sb-section {
    font-size: 0.65rem !important;
    font-weight: 700 !important;
    letter-spacing: 0.09em;
    text-transform: uppercase;
    color: var(--dim) !important;
    margin: 1.1rem 0 0.55rem !important;
    padding-bottom: 0.35rem;
    border-bottom: 1px solid var(--border);
}

[data-testid="stSidebar"] hr {
    border-color: var(--border) !important;
    margin: 0.85rem 0 !important;
}

[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {
    color: var(--text) !important;
    font-size: 0.82rem !important;
    font-weight: 600 !important;
}

[data-testid="stSidebar"] [data-baseweb="slider"] [role="slider"] {
    background: var(--accent) !important;
}

[data-testid="stSidebarNav"] {
    padding: 0.25rem 0 0.75rem;
    border-bottom: 2px solid var(--border);
    margin-bottom: 0.25rem;
}

[data-testid="stSidebarNav"] a {
    color: var(--muted) !important;
    font-size: 0.875rem !important;
    font-weight: 600 !important;
    border-radius: 8px !important;
    padding: 0.5rem 0.75rem !important;
    margin: 0.12rem 0 !important;
    border: 1px solid transparent !important;
}

[data-testid="stSidebarNav"] a:hover {
    background: #F4F4F5 !important;
    color: var(--text) !important;
    border-color: var(--border) !important;
}

[data-testid="stSidebarNav"] a[aria-current="page"] {
    background: var(--accent-soft) !important;
    color: var(--accent) !important;
    font-weight: 700 !important;
    border-color: color-mix(in srgb, var(--accent) 35%, transparent) !important;
    border-left: 3px solid var(--accent) !important;
    padding-left: calc(0.75rem - 2px) !important;
}

.sb-brand {
    padding: 0 0 1rem;
    margin-bottom: 0.25rem;
}
.sb-brand-name {
    font-size: 1.15rem;
    font-weight: 800;
    color: var(--text);
    letter-spacing: -0.03em;
}
.sb-brand-tag {
    font-size: 0.78rem;
    font-weight: 600;
    color: var(--dim);
    margin-top: 0.15rem;
}
.sb-lot-badge {
    display: inline-block;
    margin-top: 0.6rem;
    font-size: 0.68rem;
    font-weight: 700;
    letter-spacing: 0.02em;
    color: #FFFFFF;
    background: var(--lot-color, var(--accent));
    padding: 0.25rem 0.6rem;
    border-radius: 5px;
}

/* ── Composants ──────────────────────────────────────────────────────── */
[data-testid="stMetric"] {
    background: var(--surface);
    border: 2px solid var(--border);
    border-top: 3px solid var(--accent);
    border-radius: var(--radius);
    padding: 1rem 1.1rem;
    box-shadow: var(--shadow);
}
[data-testid="stMetricLabel"] {
    font-size: 0.72rem !important;
    font-weight: 700 !important;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: var(--dim) !important;
}
[data-testid="stMetricValue"] {
    font-size: 1.4rem !important;
    font-weight: 800 !important;
    color: var(--text) !important;
    letter-spacing: -0.03em;
}

.stTabs [data-baseweb="tab-list"] {
    gap: 0.2rem;
    background: var(--surface);
    border: 2px solid var(--border);
    border-radius: var(--radius);
    padding: 0.25rem;
}
.stTabs [data-baseweb="tab"] {
    border-radius: 7px;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--dim);
    padding: 0.45rem 1rem;
}
.stTabs [aria-selected="true"] {
    background: var(--accent) !important;
    color: #FFFFFF !important;
    font-weight: 700 !important;
}

.stButton > button[kind="primary"] {
    background: var(--accent) !important;
    border: 2px solid var(--accent) !important;
    color: #FFFFFF !important;
    font-weight: 700 !important;
    border-radius: 8px !important;
}
.stButton > button[kind="primary"]:hover {
    filter: brightness(1.08);
}

.stButton > button[kind="secondary"] {
    border: 2px solid var(--border) !important;
    color: var(--text) !important;
    font-weight: 600 !important;
    border-radius: 8px !important;
}

div[data-testid="stAlert"] {
    border-radius: var(--radius);
    border: 2px solid var(--border);
    font-weight: 500;
}

[data-testid="stCode"] {
    border: 2px solid var(--border);
    border-radius: 8px;
    background: #F4F4F5;
}
[data-testid="stCode"] pre {
    font-family: var(--mono) !important;
    font-size: 0.78rem !important;
    color: var(--text) !important;
}

[data-testid="stExpander"] {
    background: var(--surface);
    border: 2px solid var(--border);
    border-radius: var(--radius);
}

[data-testid="stFileUploader"] {
    background: var(--surface);
    border: 2px dashed var(--border-strong);
    border-radius: var(--radius);
}

[data-testid="stRadio"] label,
[data-testid="stCheckbox"] label {
    font-weight: 500 !important;
    color: var(--muted) !important;
}

hr { border-color: var(--border) !important; opacity: 1 !important; }

#MainMenu, footer { visibility: hidden; }
</style>
"""


def inject_styles(lot: int | None = None) -> None:
    accent = LOT_COLOR.get(lot, "#2563EB") if lot else "#2563EB"
    soft = LOT_SOFT.get(lot, "#DBEAFE") if lot else "#DBEAFE"
    lot_var = f"--lot-color: {accent};" if lot else ""
    st.markdown(
        _CSS + f"<style>:root {{ --accent: {accent}; --accent-soft: {soft}; {lot_var} }}</style>",
        unsafe_allow_html=True,
    )


def sidebar_brand(lot: int | None = None) -> None:
    badge = ""
    if lot:
        badge = f'<div class="sb-lot-badge">Lot {lot} — {LOT_LABEL[lot]}</div>'
    st.sidebar.markdown(
        f"""
        <div class="sb-brand">
            <div class="sb-brand-name">ASG-2026</div>
            <div class="sb-brand-tag">Pipeline d'assemblage de novo</div>
            {badge}
        </div>
        """,
        unsafe_allow_html=True,
    )


def sidebar_section(title: str) -> None:
    st.sidebar.markdown(f'<div class="sb-section">{title}</div>', unsafe_allow_html=True)


def page_header(title: str, subtitle: str = "") -> None:
    sub = f'<p class="sub">{subtitle}</p>' if subtitle else ""
    st.markdown(
        f'<div class="page-header"><h2>{title}</h2>{sub}</div>',
        unsafe_allow_html=True,
    )


def module_card(
    number: int,
    title: str,
    description: str,
    tags: list[str],
    lot: int | None = None,
) -> str:
    color = LOT_COLOR.get(lot, "#2563EB") if lot else "#2563EB"
    soft = LOT_SOFT.get(lot, "#DBEAFE") if lot else "#DBEAFE"
    return f"""
    <div class="dash-card" style="--card-accent:{color};--card-soft:{soft}">
        <div class="dash-card-top">
            <span class="dash-card-num">LOT {number:02d}</span>
        </div>
        <h3>{title}</h3>
        <p>{description}</p>
        <div class="dash-card-tags">{' · '.join(tags)}</div>
    </div>
    """


def stat_line(label: str, value: str) -> str:
    return f'<div class="stat-line"><span>{label}</span><b>{value}</b></div>'


def align_legend(offset: int, identity: float) -> str:
    return (
        f"Match · Mismatch · Gap — "
        f"Offset {offset} · Identité {identity * 100:.1f} %"
    )


def lot_chart(lot: int) -> dict:
    return {**CHART, "accent": LOT_COLOR[lot]}
