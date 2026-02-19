import itertools
import matplotlib.pyplot as plt
import upsetplot
from upsetplot import from_memberships
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.stats import norm, poisson

from .discretize_make_sets import (create_sets_info,
create_sets_info_twocol, create_sets_info_twocol_meta)

def create_sankey_diagram_ordered(df, vastr, possible_sets, column_prefix_dict,nbins=3):
    labellist = [f"{va}_{ps}" for va in vastr for ps in possible_sets]
    print(labellist)
    print(len(labellist), len(vastr))
    column_prefix = list(column_prefix_dict.keys())[0]
    colorlist = px.colors.qualitative.Plotly[:len(possible_sets)] * len(vastr)
    pairs = list(itertools.chain(
        *[itertools.product(range(i * len(possible_sets), (i + 1) * len(possible_sets)),
                            range((i + 1) * len(possible_sets), (i + 2) * len(possible_sets)))
          for i in range(len(vastr) - 1)]
    ))


    source = [pairs[i][0] for i in range(len(pairs))]
    target = [pairs[i][1] for i in range(len(pairs))]
    
    print("source",source)

    print("target",target)

    # print(df.columns)
    print([col for col in df.columns if column_prefix in col])
    print([np.unique(df[col].values) for col in df.columns if column_prefix in col])
    sizes = []
    possible_vals = np.arange(0, len(possible_sets))
    for va1, va2 in zip(vastr[:-1], vastr[1:]):
        for i, ps1 in enumerate(possible_vals):
            for j, ps2 in enumerate(possible_vals):

                # print(ps1,ps2,f'{va1}_{column_prefix}',f'{va2}_{column_prefix}')
                # print((df[f'{column_prefix}_{va1}'] == ps1).sum())
                # print((df[f'{column_prefix}_{va2}'] == ps2).sum())
                
                size = ((df[f'{va1}_{column_prefix}'] == ps1) & (df[f'{va2}_{column_prefix}'] == ps2)).sum()
                if size > 0:

                    sizes.append(size)
                else:
                    sizes.append(1)
    print(sizes)
    print(len(sizes))

    ch_sizes = []
    for va in vastr:
        for ps in possible_vals:
            # print(va,f'{column_prefix}_{va}',ps,(df[f'{column_prefix}_{va}'] == ps).sum())
            ch_sizes.append((df[f'{va}_{column_prefix}'] == ps).sum())

    # ch_sizes = list([((df[f'{column_prefix}_{va}'] == ps).sum() for ps in possible_sets for va in vastr )])


    print(ch_sizes)
    posxlist = np.repeat(np.linspace(0.001, 0.999, len(vastr)), len(possible_sets))
    chunk_size = len(possible_sets)
    posylist = []
    shift = 0.0
    for i in range(len(vastr)):
        chunk = ch_sizes[i * chunk_size:(i + 1) * chunk_size]
        total_size = sum(chunk)
        print("va", vastr[i], chunk, total_size)

        # Ensure the order is low, medium, high
        sorted_chunk = sorted(chunk)
        
        # Calculate the cumulative sum of the normalized sizes
        cumulative_sum = np.cumsum([(1 / (nbins + 1)) * (1 - size / total_size) for size in sorted_chunk])
        
        # posylist.extend(cumulative_sum)
        posylist.extend([0.25-shift,0.5-shift,0.75-shift])
        shift +=0.2

    print("posxlist", posxlist)
    print("posylist", posylist)

    # Rescale sizes to fit within the figure dimensions
    max_size = max(sizes)
    rescaled_sizes = [size / max_size * 100 for size in sizes]

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=[""] * len(labellist),  # Empty labels
            color=colorlist,
            x=posxlist,
            y=posylist
        ),
        link=dict(
            source=source,
            target=target,
            value=rescaled_sizes
        ),
        arrangement='snap'
    )])

    # Add va values at the bottom of each set
    for i, va in enumerate(vastr):
        fig.add_annotation(
            x=posxlist[i * len(possible_sets) + len(possible_sets) // 2],
            y=-0.1,  # Position below the plot
            text=f"v<sub>a</sub>={va[2:]}",
            showarrow=False,
            font=dict(size=20)
        )

    # Set xlim and ylim
    fig.update_xaxes(range=[0, 1], showticklabels=False)
    fig.update_yaxes(range=[-0.11, 1], showticklabels=False)
    
    # Remove xlabels, ylabels, xticks, and yticks
    fig.update_layout(
        xaxis=dict(showticklabels=False,showgrid=False),
        yaxis=dict(showticklabels=False,showgrid=False),
    )

    # Add a global color legend for low, medium, and high
    legend_labels = ["Low", "Medium", "High"]
    for i, label in enumerate(legend_labels):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=12, color=colorlist[i]),
            legendgroup=label,
            showlegend=True,
            name=label,
            textfont=dict(size=20)  # Increase the font size to 20
        ))

    fig.update_layout(
        title_text="Sankey Diagram of " + column_prefix_dict[column_prefix],
        font_size=18,
        width=800,
        height=600,
        showlegend=True,
        margin=dict(l=50, r=50, t=50, b=150)  # Adjust margins to prevent overflow
    )
    fig.show()


def generate_upset_plot(df, cva, nbins=4):
    columns_to_extract = ["Chrom1-id", "Chrom2-id", "Domain1",
                           "Domain1-points-to", "Domain2",
                           "Domain2-points-to","InteractionStrength_discretized"] + [col for col in df.columns if cva in col]
    
    # print(df.columns)
    print("cols extract",columns_to_extract)
    new_df = df[columns_to_extract]
    # print(new_df.columns)
    new_sets_df,possible_sets = create_sets_info(new_df, cva + "_MaxIntermingling_discretized",
                                                  "InteractionStrength_discretized",nbins=nbins)

    print(new_sets_df.columns)
    up_df_by_set = from_memberships(new_sets_df.sets_info.str.split(","), data=new_sets_df)
    # Rearrange the columns in the order of appearance in possible_sets
    # ordered_columns = ['chrom', 'patch', 'globindex', 'pointer', 'IAD']
    # for pair in possible_sets:
    #     for col in pair.split(","):
    #         print(col, up_df_by_set.columns)
    #         if col in new_sets_df.columns:
    #             ordered_columns.append(col)
    #             print(ordered_columns)
    # # print([col for pair in possible_sets for col in new_sets_df.columns if pair.split(",")[0] in col or pair.split(",")[1] in col])
    # up_df_by_set = up_df_by_set[ordered_columns]
    
    upsetplot.UpSet(up_df_by_set,
                    show_counts=True,
                    sort_categories_by='input',
                    # sort_by='input',
                    include_empty_subsets=False).plot()
    plt.annotate(f'$v_a$: 0{cva[3:]}', xy=(-0.85, 0.85), xycoords='axes fraction', fontsize=18, verticalalignment='top')
    plt.tight_layout()


def generate_upset_plot_w_prefix(df, annot_dict,cva, 
                                 col_prefix1,col_prefix2,nbins=4):
    columns_to_extract = ["Chrom1-id", "Chrom2-id", "Specific","InteractionStrength","InteractionStrength_discretized"] + [col for col in df.columns if cva in col]
    
    # print(df.columns)
    print("cols extract",columns_to_extract)
    new_df = df[columns_to_extract]
    # print(new_df.columns)
    new_sets_df,possible_sets = create_sets_info(new_df, cva +col_prefix1,
                                                  col_prefix2,nbins=nbins)

    print(new_sets_df.columns)
    up_df_by_set = from_memberships(new_sets_df.sets_info.str.split(","), data=new_sets_df)
    # Rearrange the columns in the order of appearance in possible_sets
    # ordered_columns = ['chrom', 'patch', 'globindex', 'pointer', 'IAD']
    # for pair in possible_sets:
    #     for col in pair.split(","):
    #         print(col, up_df_by_set.columns)
    #         if col in new_sets_df.columns:
    #             ordered_columns.append(col)
    #             print(ordered_columns)
    # # print([col for pair in possible_sets for col in new_sets_df.columns if pair.split(",")[0] in col or pair.split(",")[1] in col])
    # up_df_by_set = up_df_by_set[ordered_columns]
    
    upsetplot.UpSet(up_df_by_set,
                    show_counts=True,
                    sort_categories_by='input',
                    # sort_by='input',
                    include_empty_subsets=False).plot()
    for key, value in annot_dict.items():
        l = value[1]
        if cva[:l] == key:
            strpt1 = value[0]
            strpt2 = cva[l:]
            break
    print(strpt1 +':' +strpt2)
    plt.annotate(strpt1 +':' +strpt2, xy=(-0.85, 0.85), xycoords='axes fraction', fontsize=18, verticalalignment='top')
    plt.tight_layout()


def generate_upset_plot_w_prefix_twocol(df, cva, col_prefix1,col_prefix2,nbins1=3,nbins2=2):
    columns_to_extract = ["Chrom1-id", "Chrom2-id", "Specific","InteractionStrength","InteractionStrength_discretized"] + [col for col in df.columns if cva in col]
    
    # print(df.columns)
    print("cols extract",columns_to_extract)
    new_df = df[columns_to_extract]
    # print(new_df.columns)
    if "Specific" in col_prefix2 or "InteractionStrength" in col_prefix2:
        new_sets_df,possible_sets = create_sets_info_twocol_meta(new_df, cva +col_prefix1,
                                                  col_prefix2,nbins1=nbins1,nbins2=nbins2)
    else:
        new_sets_df,possible_sets = create_sets_info_twocol(new_df, cva +col_prefix1,
                                                  cva + col_prefix2,nbins1=nbins1,nbins2=nbins2)

    print(new_sets_df.columns)
    up_df_by_set = from_memberships(new_sets_df.sets_info.str.split(","), data=new_sets_df)
    # Rearrange the columns in the order of appearance in possible_sets
    # ordered_columns = ['chrom', 'patch', 'globindex', 'pointer', 'IAD']
    # for pair in possible_sets:
    #     for col in pair.split(","):
    #         print(col, up_df_by_set.columns)
    #         if col in new_sets_df.columns:
    #             ordered_columns.append(col)
    #             print(ordered_columns)
    # # print([col for pair in possible_sets for col in new_sets_df.columns if pair.split(",")[0] in col or pair.split(",")[1] in col])
    # up_df_by_set = up_df_by_set[ordered_columns]
    
    upsetplot.UpSet(up_df_by_set,
                    show_counts=True,
                    sort_categories_by='input',
                    # sort_by='input',
                    include_empty_subsets=False).plot()
    plt.annotate(f'$v_a$: 0{cva[3:]}', xy=(-0.85, 0.85), xycoords='axes fraction', fontsize=18, verticalalignment='top')
    plt.tight_layout()