streamlit as st
import pandas as pd
import requests
from bs4 import BeautifulSoup

# --- Expression Dictionary ---
default_expression = {
    "CD3": "high", "CD4": "high", "CD8": "high", "CD19": "high",
    "CD56": "medium", "CD14": "high", "HLA-DR": "medium", "CD123": "medium",
    "CD11c": "low", "CD25": "medium", "FoxP3": "low", "PD-1": "low"
}

# --- Scraper ---
def fetch_markers_cellmarker(cell_type):
    url = f"https://www.cellmarker.com/?s={cell_type.replace(' ', '+')}"
    r = requests.get(url, timeout=10)
    soup = BeautifulSoup(r.text, "html.parser")
    markers = []
    for h2 in soup.find_all('h2'):
        if 'Marker' in h2.text or 'Markers' in h2.text:
            ul = h2.find_next_sibling('ul')
            if ul:
                for li in ul.find_all('li'):
                    parts = li.text.split("–")
                    if len(parts) == 2:
                        marker, desc = parts[0].strip(), parts[1].strip()
                        if marker.startswith('CD'):
                            markers.append({"Marker": marker, "Description": desc})
            break
    return markers

def annotate_expression(markers):
    for m in markers:
        m["Expression"] = default_expression.get(m["Marker"], "medium")
    return markers

def assign_fluorochromes(cell_name, markers):
    fluorochromes = {
        "PE": 10, "APC": 9, "BV421": 8, "BUV737": 8,
        "BV510": 7, "BUV395": 7, "FITC": 4, "PerCP-Cy5.5": 3,
        "APC-Cy7": 6, "PE-Cy7": 6, "AF700": 5, "BV605": 6,
        "BV650": 6, "BV786": 5, "BB515": 4, "Zombie Aqua": 10
    }
    used = []
    panel = []
    for m in markers:
        best_score = float("inf")
        best_fluor = None
        for f, b in fluorochromes.items():
            if f in used:
                continue
            score = (10 - b) if m["Expression"] == "low" else (b if m["Expression"] == "high" else abs(6 - b))
            if score < best_score:
                best_score = score
                best_fluor = f
        used.append(best_fluor)
        panel.append({"Marker": m["Marker"], "Expression": m["Expression"], "Fluorochrome": best_fluor})
    return pd.DataFrame(panel)

def generate_gating_strategy(markers):
    return [
        "1. Live/Dead gating to exclude dead cells.",
        "2. FSC vs SSC to gate lymphocytes or monocytes.",
        "3. CD45+ gating for leukocyte lineage.",
        "4. Lineage gating (e.g., CD3, CD4/CD8, CD19).",
        "5. Subset-specific gates (activation markers, cytokines, etc.)."
    ]

def generate_fmo_plan(markers):
    return [m["Marker"] for m in markers if m["Expression"] in ["low", "medium"]]

def generate_compensation_guidance(markers):
    return "Use single-stained controls for each fluorochrome to generate compensation matrix."

# --- Streamlit App ---
st.title("🧬 Auto Flow Cytometry Panel Builder")
cell_type = st.text_input("Enter human immune cell type (e.g. 'T cell', 'NK cell')")

if st.button("Design Panel"):
    with st.spinner("Fetching markers and designing your panel..."):
        raw_markers = fetch_markers_cellmarker(cell_type)
        if not raw_markers:
            st.error("No markers found for that cell type.")
        else:
            annotated = annotate_expression(raw_markers)
            panel_df = assign_fluorochromes(cell_type, annotated)
            st.subheader("🔬 Fluorochrome Assignment")
            st.dataframe(panel_df)

            st.subheader("🧭 Gating Strategy")
            for step in generate_gating_strategy(annotated):
                st.write("- " + step)

            st.subheader("🔍 FMO Control Suggestions")
            st.write(", ".join(generate_fmo_plan(annotated)))

            st.subheader("🎛️ Compensation Guidance")
            st.write(generate_compensation_guidance(annotated))

            st.download_button("⬇️ Download CSV", panel_df.to_csv(index=False), "panel.csv")
