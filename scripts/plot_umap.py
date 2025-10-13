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
st.title('Building Block Cluster Explorer')

# initialize session state variables
if 'plot_version' not in st.session_state:
    st.session_state.plot_version = 0

if 'sampled_indices' not in st.session_state:
    st.session_state.sampled_indices = []

if 'compound_index' not in st.session_state:
    st.session_state.compound_index = 0

# load data
acids_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_acids.parquet')
amines_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_amines.parquet')

# add sidebar
st.sidebar.header('Settings')

# add building block type toggle buttons
bb_type = st.sidebar.radio(
    'Functional group:',
    ['Amines', 'Acids']
)

# reset sampled indices if BB type changes
# Clear sampled indices if BB type changes
if 'previous_bb_type' not in st.session_state:
    st.session_state.previous_bb_type = bb_type
elif st.session_state.previous_bb_type != bb_type:
    st.session_state.sampled_indices = []
    st.session_state.previous_bb_type = bb_type

# select appropriate bb df from toggle
if bb_type == 'Amines':
    df = amines_df.set_index('eMolecules ID')
else:
    df = acids_df.set_index('eMolecules ID')

# add sampling inputs and controls
st.sidebar.markdown('---')
st.sidebar.subheader('Building Block Sampler')
sample_size = st.sidebar.number_input(
    'Sample size:',
    min_value=1,
    max_value=len(df),
    value=100,
    step=1
)

#col_sample1, col_sample2 = st.sidebar.columns(2)
#with col_sample1:
if st.sidebar.button('📊 Sample Points', width='stretch'):
    # Randomly sample points
    #sampled_indices = df.sample(n=min(sample_size, len(df))).index.tolist()
    #st.session_state.sampled_indices = sampled_indices
    #st.rerun()

    # get cluster proportions
    cluster_counts = df['cluster'].value_counts()
    total_samples = min(sample_size, len(df))

    sampled_dfs = []

    for cluster in cluster_counts.index:
        cluster_df = df[df['cluster'] == cluster]

        # compute samples per cluster (ensure at least 1 per cluster)
        cluster_proportion = len(cluster_df)/len(df)
        n_samples = max(1, round(cluster_proportion*total_samples))

        # ensure n_samples is valid
        n_samples = min(n_samples, len(cluster_df))

        # sample from cluster
        sampled_cluster = cluster_df.sample(n=n_samples)
        sampled_dfs.append(sampled_cluster)

    # combine all sampled dfs
    combined_samples = pd.concat(sampled_dfs)

    # ensure total num samples is valid (randomly drop if too many)
    if len(combined_samples) > total_samples:
        combined_samples = combined_samples.sample(n=total_samples)
    
    st.session_state.sampled_indices = combined_samples.index.tolist()
    st.rerun()

#with col_sample2:
if st.sidebar.button('Clear Sample', width='stretch'):
    st.session_state.sampled_indices = []
    st.rerun()

# Show sampling info
if 'sampled_indices' in st.session_state and st.session_state.sampled_indices:
    st.sidebar.success(f'✓ {len(st.session_state.sampled_indices)} points sampled')

# display summary
st.subheader('Dataset summary')
col1, col2, col3 = st.columns(3)
with col1:
    st.metric('Total compounds', len(df))
with col2:
    st.metric('Number of clusters', df['cluster'].nunique())
with col3:
    st.metric('Selected functional group', bb_type)

# show cluster stats
st.markdown('**Cluster details**')
cluster_counts = df['cluster'].value_counts().sort_index()
cluster_counts = cluster_counts.sort_index()
cluster_df = pd.DataFrame({
    'Cluster': cluster_counts.index,
    'Count': cluster_counts.values,
    'Percentage (%)': (cluster_counts.values/len(df)*100).round(1)
})
st.dataframe(cluster_df, hide_index=True, width='stretch')

st.divider()

# convert cluster to categorical type
df['cluster'] = df['cluster'].astype(str)

# format price column
df['Price_usd'] = df['Price'].apply(lambda x: f'${x:.2f}')

# set color palette
palette = q.Plotly

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
    },
    title=f'UMAP projection for {bb_type} building blocks',
    color_discrete_sequence=palette,
    height=600
)

# update style and markers
fig.update_traces(
    marker=dict(
        size=10,
        line=dict(width=1, color='white')
    ),
)

# update layout
fig.update_layout(
    plot_bgcolor='white',
    paper_bgcolor='white',
    font=dict(size=12),
    xaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray'),
    yaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray')
)

# Add highlighted layer for sampled points if any exist
if st.session_state.sampled_indices:
    # filter to only indices that exist in current df
    valid_indices = [idx for idx in st.session_state.sampled_indices if idx in df.index]

    #sampled_df = df.loc[st.session_state.sampled_indices]
    if valid_indices:
        sampled_df = df.loc[valid_indices]

        # Add sampled points as a separate trace with highlighting
        fig.add_scatter(
            x=sampled_df['UMAP1'],
            y=sampled_df['UMAP2'],
            mode='markers',
            marker=dict(
                size=15,
                color='DarkSlateGrey',
                #line=dict(width=3, color='red')  # Bold red outline
            ),
            name='Sampled',
            showlegend=True,
            hoverinfo='skip'  # Don't show hover for this layer
        )

        #### z order #####


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
    if st.button('🔄  Reset Plot Selection'):
        # Increment plot version to force a fresh render with all traces visible
        st.session_state.plot_version += 1
        st.session_state.compound_index = 0
        st.rerun()    

with col_details:
    # Use the dynamic key to access the selection
    plot_key = f'scatter_plot_{st.session_state.plot_version}'
    
    # Display molecule card when point is clicked
    if plot_key in st.session_state and st.session_state[plot_key]:
        selection = st.session_state[plot_key].get('selection')
        if selection and 'points' in selection and len(selection['points']) > 0:
            
            num_selected = len(selection['points'])
            
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

# Show sampled compounds details
with st.expander(f'Sampled Building Block Details ({len(st.session_state.sampled_indices)})'):
    if st.session_state.sampled_indices:
        #st.subheader(f'Sampled Compounds ({len(st.session_state.sampled_indices)})')
        
        sampled_df = df.loc[st.session_state.sampled_indices]
        st.dataframe(sampled_df[[
            'eMolecules SKU',
            'SMILES',
            'Price',
            'Packsize',
            'Supplier Name',
            'Tier',
            'cluster'
        ]], use_container_width=True)
        
        # Add download button for sampled data
        csv = sampled_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Sampled Data as CSV",
            data=csv,
            file_name=f'sampled_{bb_type.lower()}_{len(st.session_state.sampled_indices)}_compounds.csv',
            mime='text/csv',
        )