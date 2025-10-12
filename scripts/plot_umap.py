import streamlit as st
import pandas as pd
import plotly.express as px
from plotly.colors import qualitative as q
from streamlit_plotly_events import plotly_events
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# helper function to render chemical structure
def render_molecule(smiles):
    """Render a molecule from SMILES string as SVG"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg
    except:
        return None

# configure page
st.set_page_config(page_title='UMAP Cluster Explorer', layout='wide')
st.title('UMAP Cluster Visualization')
st.markdown('Interactive 2D visualization of building block clusters')

# load data
acids_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_acids.parquet')
amines_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_amines.parquet')

# add sidebar for BB toggle
st.sidebar.header('Settings')
bb_type = st.sidebar.radio(
    'Building block type:',
    ['Amines', 'Acids']
)

# select appropriate bb df
if bb_type == 'Amines':
    df = amines_df
    #color_scheme = 'Blues'
else:
    df = acids_df
    #color_scheme = 'Reds'

# display stats
col1, col2, col3 = st.columns(3)
with col1:
    st.metric('Total compounds', len(df))
with col2:
    st.metric('Number of clusters', df['cluster'].nunique())
with col3:
    st.metric('Selected type', bb_type)

# convert cluster to categorical type
df['cluster'] = df['cluster'].astype(str)

# format price column
df['Price_usd'] = df['Price'].apply(lambda x: f'${x:.2f}')

# set color palette
palette = q.Plotly  # other options: q.Set3, q.Bold, q.Dark24, q.Light24, q.Vivid

# Initialize plot version in session state (used to force plot refresh)
if 'plot_version' not in st.session_state:
    st.session_state.plot_version = 0

# create UMAP plot
fig = px.scatter(
    df,
    x='UMAP1',
    y='UMAP2',
    color='cluster',
    hover_data={
        'UMAP1': False,
        'UMAP2': False,
        'eMolecules SKU': True,
        'SMILES': True,
        'Price_usd': True,
        'Packsize': True,
        'cluster': True
    },
    title=f'UMAP projection with clustering for {bb_type} building blocks',
    color_discrete_sequence=palette,
    height=600
)

# update style and markers
fig.update_traces(
    marker=dict(
        size=12,
        line=dict(width=1, color='white')
    ),
    hovertemplate='<br>'.join([
        'eMolecules SKU: %{customdata[0]}',
        'SMILES: %{customdata[1]}',
        'Price: %{customdata[2]}',
        'Packsize: %{customdata[3]}',
        'Cluster: %{customdata[4]}',
        '<extra></extra>'
    ])
)

# update layout
fig.update_layout(
    plot_bgcolor='white',
    paper_bgcolor='white',
    font=dict(size=12),
    xaxis=dict(showgrid=True, gridwidth=1,gridcolor='lightgray'),
    yaxis=dict(showgrid=True, gridwidth=1,gridcolor='lightgray')
)

# Create two columns for plot and details panel
col_plot, col_details = st.columns([2, 1])

with col_plot:
    # display plot with a dynamic key to force refresh when needed
    selected = st.plotly_chart(
        fig, 
        use_container_width=True, 
        key=f'scatter_plot_{st.session_state.plot_version}',  # Dynamic key
        on_select='rerun'
    )

    # Add reset button above the plot
    if st.button('🔄  Reset Plot'):
        # Increment plot version to force a fresh render with all traces visible
        st.session_state.plot_version += 1
        st.session_state.compound_index = 0
        st.rerun()

    st.markdown('---')
    st.markdown('Double click on a cluster label to view single cluster. ' \
    'Single click to view or hide additional clusters')

with col_details:
    # Use the dynamic key to access the selection
    plot_key = f'scatter_plot_{st.session_state.plot_version}'
    
    # Display molecule card when point is clicked
    if plot_key in st.session_state and st.session_state[plot_key]:
        selection = st.session_state[plot_key].get('selection')
        if selection and 'points' in selection and len(selection['points']) > 0:
            
            num_selected = len(selection['points'])
            
            # Initialize the current index in session state
            if 'compound_index' not in st.session_state:
                st.session_state.compound_index = 0
            
            # Reset index if it's out of bounds (happens when switching selections)
            if st.session_state.compound_index >= num_selected:
                st.session_state.compound_index = 0
            
            # Navigation controls for multiple selections
            if num_selected > 1:
                st.subheader(f'Selected Compounds ({num_selected})')
                
                col_nav1, col_nav2, col_nav3 = st.columns([1, 2, 1])
                
                with col_nav1:
                    if st.button('← Prev', use_container_width=True):
                        st.session_state.compound_index = (st.session_state.compound_index - 1) % num_selected
                        st.rerun()
                
                with col_nav2:
                    st.markdown(f"<h4 style='text-align: center;'>Compound {st.session_state.compound_index + 1} of {num_selected}</h4>", unsafe_allow_html=True)
                
                with col_nav3:
                    if st.button('Next →', use_container_width=True):
                        st.session_state.compound_index = (st.session_state.compound_index + 1) % num_selected
                        st.rerun()
                
                st.markdown('---')
            else:
                st.subheader('Selected Compound')
            
            # Get the current compound to display
            point_info = selection['points'][st.session_state.compound_index]
            emolecules_sku = point_info['customdata'][0]
            selected_row = df[df['eMolecules SKU'] == emolecules_sku].iloc[0]
            
            # Chemical structure
            st.markdown('**Chemical Structure**')
            smiles = selected_row['SMILES']
            svg = render_molecule(smiles)
            
            if svg:
                st.image(svg, use_container_width=True)
            else:
                st.error('Unable to render molecule structure')
            
            st.markdown('---')
            # Compound details
            st.markdown('**Compound Details**')
            
            details = {
                'eMolecules SKU': selected_row['eMolecules SKU'],
                'eMolecules ID': selected_row['eMolecules ID'],
                'SMILES': selected_row['SMILES'],
                'CAS': selected_row.get('CAS', 'N/A'),
                'Cluster': selected_row['cluster'],
                'Price': f"${selected_row['Price']:.2f}" if pd.notna(selected_row['Price']) else 'N/A',
                'Packsize': selected_row['Packsize'],
                'Tier': selected_row.get('Tier', 'N/A'),
                'Supplier Name': selected_row.get('Supplier Name', 'N/A'),
                'Supplier Catalog Number': selected_row.get('Supplier Catalog Number', 'N/A')
            }
            
            for key, value in details.items():
                st.text(f"{key}: {value}")

        else:
            # Reset index when nothing is selected
            st.session_state.compound_index = 0
            st.info('👈 Click on a point in the plot to view compound details')
    else:
        # Show placeholder when nothing selected
        st.session_state.compound_index = 0
        st.info('👈 Click on a point in the plot to view compound details')

st.divider()

# show sample data
with st.expander('Catalog Details'):
    st.dataframe(df[[
        'eMolecules ID',
        'eMolecules SKU',
        'SMILES',
        'CAS',
        'Tier',
        'Price',
        'Packsize',
        'Supplier Name',
        'Supplier Catalog Number'
    ]].head(100), width='stretch')

# show cluster stats
st.subheader('Cluster details')
cluster_counts = df['cluster'].value_counts().sort_index()
cluster_df = pd.DataFrame({
    'Cluster': cluster_counts.index,
    'Count': cluster_counts.values,
    'Percentage': (cluster_counts.values/len(df)*100).round(1)
})
st.dataframe(cluster_df, hide_index=True, width='stretch')