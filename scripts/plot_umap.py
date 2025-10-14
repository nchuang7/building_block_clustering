import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import re
from plotly.colors import qualitative as q
from streamlit_plotly_events import plotly_events
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

##### CONFIGURATION #####

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

if 'stratified_indices' not in st.session_state:
    st.session_state.stratified_indices = []

if 'manual_indices' not in st.session_state:
    st.session_state.manual_indices = []

if 'compound_index' not in st.session_state:
    st.session_state.compound_index = 0

if 'highlight_skus' not in st.session_state:
    st.session_state.highlight_skus = []

if 'last_selection_skus' not in st.session_state:
    st.session_state.last_selection_skus = []

# load data
acids_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_acids.parquet')
amines_df = pd.read_parquet('/Users/nataliechuang/Documents/Personal Projects/Satomic/building_block_clustering/data/sampled_amines.parquet')

##### SIDEBAR #####

# add sidebar
st.sidebar.header('Settings')

# add building block type toggle buttons
bb_type = st.sidebar.radio(
    'Functional group:',
    ['Amines', 'Acids']
)

# reset sampled indices if BB type changes
if 'previous_bb_type' not in st.session_state:
    st.session_state.previous_bb_type = bb_type
elif st.session_state.previous_bb_type != bb_type:
    st.session_state.sampled_indices = []
    st.session_state.stratified_indices = []
    st.session_state.manual_indices = []
    st.session_state.previous_bb_type = bb_type

# toggle appropriate bb
if bb_type == 'Amines':
    df = amines_df.set_index('eMolecules ID')
else:
    df = acids_df.set_index('eMolecules ID')

# add sampling inputs and controls
st.sidebar.markdown('---')
st.sidebar.subheader('Building Block Sampler')

# update sampled indices
st.session_state.sampled_indices = list(set(st.session_state.stratified_indices + st.session_state.manual_indices))

## STRATIFIED SAMPLING ##

# sample size input
st.sidebar.markdown('**Stratified Sampling**')
sample_size = st.sidebar.number_input(
    'Select sample size:',
    min_value=1,
    max_value=len(df),
    value=100,
    step=1,
    key='sample_size_input'
)

# add sample button
if st.sidebar.button('📊 Add Sample', width='stretch', key='add_sample'):
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
        sampled_cluster = cluster_df.sample(n=n_samples, random_state=7)
        sampled_dfs.append(sampled_cluster)

    # combine all sampled dfs
    combined_samples = pd.concat(sampled_dfs)

    # ensure total num samples is valid (randomly drop if too many)
    if len(combined_samples) > total_samples:
        combined_samples = combined_samples.sample(n=total_samples, random_state=7)

    # combine with any existing stratified samples
    existing_stratified = set(st.session_state.stratified_indices)
    new_indices = set(combined_samples.index.tolist())
    st.session_state.stratified_indices = list(existing_stratified.union(new_indices))

    # update combined sample indices
    st.session_state.sampled_indices = list(set(st.session_state.stratified_indices + st.session_state.manual_indices))
    st.rerun()

# clear sample button
if st.sidebar.button('Clear Sample', width='stretch', key='clear_sample'):
    # clear only stratified samples
    st.session_state.stratified_indices = []
    
    # update combined sample indices
    st.session_state.sampled_indices = list(set(st.session_state.manual_indices))
    st.rerun()

st.sidebar.markdown('---')

## MANUAL SAMPLING ##

# manual sample input box
st.sidebar.markdown('**Manual ID Input**')
st.sidebar.caption('Enter eMolecules IDs to *add to* sample list (comma, space, or newline separated):')
id_input = st.sidebar.text_area(
    'eMolecules IDs:',
    height=120,
    placeholder='e.g. 12345, 67890 11111\n or one per line',
    label_visibility='collapsed',
    key='id_input'
)

# add ID button
if st.sidebar.button('🎯 Add IDs', width='stretch', key='add_ids'):
    if id_input.strip():
        # parse input to handle comma, space, or newline separated
        id_strings = re.split(r'[,\s\n]+', id_input.strip())
            
        try:
            # Convert to integers
            manual_ids = [int(id_str.strip()) for id_str in id_strings if id_str.strip()]
                
            # Find matching indices in dataframe
            matching_indices = df[df.index.isin(manual_ids)].index.tolist()
                
            if matching_indices:
                # Add to manual indices (using set to avoid duplicates)
                existing_manual = set(st.session_state.manual_indices)
                new_indices = set(matching_indices)
                st.session_state.manual_indices = list(existing_manual.union(new_indices))
                
                # Update combined sampled_indices
                st.session_state.sampled_indices = list(set(st.session_state.stratified_indices + st.session_state.manual_indices))

                not_found = set(manual_ids) - set(matching_indices)
                if not_found:
                    st.sidebar.warning(f'Found {len(matching_indices)} IDs. Could not find: {len(not_found)} IDs')
                st.rerun()
            else:
                st.sidebar.error('No matching eMolecules IDs found in current dataset')
        except ValueError:
            st.sidebar.error('Invalid input. Please enter valid numeric IDs.')
    else:
       st.sidebar.warning('Please enter at least one eMolecules ID')

# clear ID button
if st.sidebar.button('Clear IDs', width='stretch', key='clear_ids'):
    # clear only manual IDs
    st.session_state.manual_indices = []

    # update combined sample indices
    st.session_state.sampled_indices = list(set(st.session_state.stratified_indices))
    st.rerun()

st.sidebar.markdown('---')

## MANUAL SAMPLE REMOVAL ##

# manual sample removal input box
st.sidebar.markdown('**Manual ID Removal**')
st.sidebar.caption('Enter eMolecules IDs to *remove from* sample list (comma, space, or newline separated):')
remove_input = st.sidebar.text_area(
    'IDs to remove:',
    height=100,
    placeholder='e.g. 12345, 67890',
    label_visibility='collapsed',
    key='remove_input'
)

# remove ID button
if st.sidebar.button('❌ Remove IDs', width='stretch', key='remove_ids'):
    if remove_input.strip():
        # Parse input
        id_strings = re.split(r'[,\s\n]+', remove_input.strip())
        
        try:
            # Convert to integers
            ids_to_remove = [int(id_str.strip()) for id_str in id_strings if id_str.strip()]
            
            # Remove from both stratified and manual lists
            st.session_state.stratified_indices = [idx for idx in st.session_state.stratified_indices if idx not in ids_to_remove]
            st.session_state.manual_indices = [idx for idx in st.session_state.manual_indices if idx not in ids_to_remove]
            
            # Update combined sampled_indices
            st.session_state.sampled_indices = list(set(st.session_state.stratified_indices + st.session_state.manual_indices))
            
            # Show feedback
            removed_count = len([idx for idx in ids_to_remove if idx not in st.session_state.sampled_indices])
            if removed_count > 0:
                st.sidebar.success(f'Removed {removed_count} IDs')
            else:
                st.sidebar.warning('None of the specified IDs were in the sample')
            
            st.rerun()
        except ValueError:
            st.sidebar.error('Invalid input. Please enter valid numeric IDs.')
    else:
        st.sidebar.warning('Please enter at least one eMolecules ID to remove')

st.sidebar.markdown('---')

# clear all sample button
if st.sidebar.button('🗑️ Clear All Samples', width='stretch', type='primary'):
    st.session_state.sampled_indices = []
    st.session_state.stratified_indices = []
    st.session_state.manual_indices = []
    st.rerun()

# Show sampling info
if 'sampled_indices' in st.session_state and st.session_state.sampled_indices:
    valid_indices = [idx for idx in st.session_state.sampled_indices if idx in df.index]
    st.sidebar.success(f'✓ {len(valid_indices)} points sampled')

##### SUMMARY INFO #####

# display dataset summary
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
st.dataframe(cluster_df, hide_index=True, width='content')

st.divider()

# dataframe formatting

# convert cluster to categorical type
df['cluster'] = df['cluster'].astype(str)

# format price column
df['Price_usd'] = df['Price'].apply(lambda x: f'${x:.2f}')

##### PLOTTING #####

## PLOT CONFIGURATION ##

# set color palette
palette = q.Plotly

# define initial plot key
plot_key = f'scatter_plot_{st.session_state.plot_version}'

# handle selection first (available from the *previous* render)
selection = None
if plot_key in st.session_state and st.session_state[plot_key]:
    selection = st.session_state[plot_key].get('selection')

    # If we have a NEW selection from the plot, update highlight_skus
    if selection and 'points' in selection and len(selection['points']) > 0:
        selected_skus = [pt['customdata'][0] for pt in selection['points']]
        st.session_state.highlight_skus = selected_skus
        # Reset compound index when NEW selection is made
        if st.session_state.get('last_selection_skus') != selected_skus:
            st.session_state.compound_index = 0
            st.session_state.last_selection_skus = selected_skus

# update highlight_skus + compound_index + current selected row
current_row = None

if st.session_state.highlight_skus:
    num_selected = len(st.session_state.highlight_skus)
    if st.session_state.compound_index >= num_selected:
        st.session_state.compound_index = 0

    # get active SKU for the detail panel
    active_sku = st.session_state.highlight_skus[st.session_state.compound_index]
    if active_sku in df['eMolecules SKU'].values:
        current_row = df[df['eMolecules SKU'] == active_sku].iloc[0]
else:
    # no selection -> clear point halo
    st.session_state.compound_index = 0

## FIGURE DEFINITION AND FORMATTING ##

# create figure for scatter plot
fig = go.Figure()

# manually add each cluster as a trace - use Scattergl for performance
for cluster in sorted(df['cluster'].unique()):
    cluster_df = df[df['cluster'] == cluster]
    fig.add_trace(go.Scattergl(
        x=cluster_df['UMAP1'],
        y=cluster_df['UMAP2'],
        mode='markers',
        name=str(cluster),
        marker=dict(
            size=10,
            color=palette[int(cluster) % len(palette)],
            line=dict(width=1, color='white')
        ),
        customdata=cluster_df[['eMolecules SKU', 'SMILES']].values,
        hovertemplate='<br>'.join([
            'eMolecules SKU=%{customdata[0]}',
            'SMILES=%{customdata[1]}',
            '<extra></extra>'
        ])
    ))

# add sampled points as last trace
if st.session_state.sampled_indices:
    valid_indices = [idx for idx in st.session_state.sampled_indices if idx in df.index]
    if valid_indices:
        sampled_df = df.loc[valid_indices]
        fig.add_trace(go.Scattergl(
            x=sampled_df['UMAP1'],
            y=sampled_df['UMAP2'],
            mode='markers',
            name='Sampled',
            marker=dict(
                size=14,
                color='black',
            ),
            customdata=sampled_df[['eMolecules SKU', 'SMILES']].values,
            hoverinfo='skip'
        ))

# update layout
fig.update_layout(
    plot_bgcolor='white',
    paper_bgcolor='white',
    font=dict(size=12),
    xaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray', title='UMAP1'),
    yaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray', title='UMAP2'),
    height=600
)

# highlight selected points on plot
if st.session_state.highlight_skus:
    highlight_df = df[df['eMolecules SKU'].isin(st.session_state.highlight_skus)]
    if not highlight_df.empty:
        fig.add_trace(go.Scattergl(
            x=highlight_df['UMAP1'],
            y=highlight_df['UMAP2'],
            mode='markers',
            name='Selected',
            hoverinfo='skip',
            marker=dict(
                size=20,
                color='rgba(0,0,0,0)',
                line=dict(width=3, color='black')
            ),
            customdata=highlight_df[['eMolecules SKU', 'SMILES']].values,
            showlegend=False
        ))

## PLOT DISPLAY AND DETAIL PANEL ##

# create two columns for plot and details panel
col_plot, col_details = st.columns([2, 1])

# display plot
with col_plot:
    st.markdown(f'### UMAP projection for {bb_type} building blocks')
    selected = st.plotly_chart(
        fig,
        use_container_width=True,
        key=plot_key,
        on_select='rerun'
    )

    # add reset button
    if st.button('🔄  Reset Plot Selection'):
        st.session_state.plot_version += 1
        st.session_state.compound_index = 0
        st.session_state.highlight_skus = []
        st.session_state.last_selection_skus = []
        st.rerun()     

# add building block structure and details
with col_details:
    if current_row is not None:
        selected_count = len(st.session_state.highlight_skus)
        if selected_count > 1:
            col_nav1, col_nav2, col_nav3 = st.columns([1,2,1])
            with col_nav1:
                if st.button('← Prev', use_container_width=True):
                    st.session_state.compound_index = (st.session_state.compound_index - 1) % selected_count
                    st.rerun()
            with col_nav2:
                st.markdown(
                    f"<h4 style='text-align: center;'>Compound {st.session_state.compound_index + 1} of {selected_count}</h4>",
                    unsafe_allow_html=True
                )
            with col_nav3:
                if st.button('Next →', use_container_width=True):
                    st.session_state.compound_index = (st.session_state.compound_index + 1) % selected_count
                    st.rerun()

        # structure display
        st.markdown('**Structure**')
        smiles = current_row['SMILES']
        svg = render_molecule(smiles)
        if svg:
            st.image(svg)
        else:
            st.error('Unable to render molecule structure')

        # detail display
        st.markdown('---')
        st.markdown('**Details**')
        details = {
            'eMolecules SKU': current_row['eMolecules SKU'],
            'SMILES': current_row['SMILES'],
            'CAS': current_row.get('CAS', 'N/A'),
            'Cluster': current_row['cluster'],
            'Price': f"${current_row['Price']:.2f}" if pd.notna(current_row['Price']) else 'N/A',
            'Packsize': current_row.get('Packsize', 'N/A'),
            'Tier': current_row.get('Tier', 'N/A'),
            'Supplier Name': current_row.get('Supplier Name', 'N/A'),
            'Supplier Catalog Number': current_row.get('Supplier Catalog Number', 'N/A')
        }
        for k, v in details.items():
            st.text(f"{k}: {v}")
    else:
        st.info('👈 Click on a point in the plot to view compound details')

st.divider()

##### SAMPLE DETAIL TABLE #####

if st.session_state.sampled_indices:
    # filter to valid indices that exist in current dataframe
    valid_sampled_indices = [idx for idx in st.session_state.sampled_indices if idx in df.index]
        
    if valid_sampled_indices:
        sampled_display_df = df.loc[valid_sampled_indices]
        sampled_display_df = sampled_display_df.rename(columns={'Price': 'Price (USD$)', 'cluster': 'Cluster'})

        st.subheader(f'Sampled Building Block Details ({len(valid_sampled_indices)})')
        st.dataframe(sampled_display_df[[
            'eMolecules SKU',
            'SMILES',
            'Price (USD$)',
            'Packsize',
            'Supplier Name',
            'Tier',
            'Cluster'
        ]], width='stretch')
        
        # Add download button for sampled data
        csv = sampled_display_df.to_csv()
        st.download_button(
            label="📥 Download Sampled Data as CSV",
            data=csv,
            file_name=f'sampled_{bb_type.lower()}_{len(valid_sampled_indices)}_building_blocks.csv',
            mime='text/csv',
        )