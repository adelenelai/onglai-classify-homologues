import mols2grid
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
import rdkit.Chem.Draw as Draw
import pandas as pd

st.title("Explore Classified Homologous Series")

st.subheader("Introduction")
st.markdown(
    "Homologous series are groups of related chemical compounds that share a common core substructure and a repeating chemical subunit. Compounds in a homologous series can be represented by a general formula or Markush structure. A basic example is alkanes (CnH2n+1): methane, ethane, propane *etc.*"
)
st.markdown(
    "Many areas of chemistry have homologous series, such as natural products and environmental chemicals (*e.g.,* High Production Volume chemicals like surfactants). Because of their structural similartiy, compounds within the same homologous series typically share similar or predictable properties, and may produce distinct analytical signatures."
)

st.subheader("Classifiying Homologous Series")
st.markdown(
    "Within most compound datasets or databases, homologous series are typically unclassified. Therefore, an algorithm was developed and implemented using the RDKit to classify homologous series."
)

st.subheader("Visualise Homologous Series Classified")
file_in = st.file_uploader("Upload CSV results:", type="csv")

#### summary stats - figure of distribution

# df = pd.read_csv(file_in)
df = pd.read_csv(
    "classified_series.csv"
)  # default test1_23.csv as trial data within Docker image?


# ##Dropdown box
# series_select = st.selectbox(
#     label='Select a classified homologous series by series_no:',
#     options=df.series_no.unique(),
#     index=0)

## Number input
series_select = st.number_input(
    label="Select a classified homologous series by series_no:",
    min_value=-3,
    max_value=max(df.series_no),
    value=0,
    step=1,
)

## Subset df depending on series_no
df_result = df[df["series_no"] == series_select]


## Get core
core_smiles = df.loc[df["series_no"] == series_select, "CanoSmiles_FinalCores"].iloc[1]
core_result = Chem.MolFromSmiles(core_smiles)


st.markdown("**Core detected**")
st.image(Draw.MolToImage(core_result), caption=core_smiles, width=250)


st.markdown("**Homologous series**")
raw_html_series = mols2grid.display(
    df_result,
    # mapping={"smiles": "SMILES", "Name": "cpd_name"},
    tooltip=["molecular_formula"],
    subset=["img", "cpd_name"],
    size=(170, 140),
    n_cols=5,
    n_rows=3,
    sort_by="molecular_formula",
)._repr_html_()


components.html(raw_html_series, height=300, width=900, scrolling=True)

st.write(df_result.loc[:, ["cpd_name", "SMILES", "InChIKey", "monoisotopic_mass"]],)

st.subheader("**Links & Citation**")
st.markdown(
    "1. ```Lai, A., Schaub, J., Steinbeck, C., Schymanski, E. L. (2022). An Algorithm to Classify Homologous Series.```*in prep*"
)
st.markdown(
    "2. [GitHub Repository](https://github.com/adelenelai/classify_homologues) containing HS classification algorithm code."
)
