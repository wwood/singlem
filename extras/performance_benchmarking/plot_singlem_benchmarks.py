"""Plot wall time and peak RAM of SingleM pipe/renew/genome across versions."""
import argparse, glob, os
import polars as pl
from plotnine import *

VERSIONS = ['0-19-0', '0-20-3', '0-21-4', 'dev-2026-04', 'dev']
LABELS = {'0-19-0': '0.19.0', '0-20-3': '0.20.3', '0-21-4': '0.21.4', 'dev-2026-04': 'dev (Apr 2026)', 'dev': 'dev (now)'}
COLORS = ['#2a78d6', '#eb6834', '#6250d6', '#1baf7a', '#eda100']
STAGES = {'singlem-no_taxonomy': 'pipe (no taxonomy)', 'singlem-renew': 'renew', 'singlem-genome': 'pipe --genome-fasta'}
SAMPLES = {'DRR083195': 'DRR083195\nIllumina', 'SRR1926147': 'SRR1926147\nIllumina',
           'SRR9224014': 'SRR9224014\nIllumina', 'ERR11563597': 'ERR11563597\nnanopore',
           'methanobacteria_genomes': 'Methanobacteria\ngenomes'}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir', default='snakedata_r232',
                    help='snakemake output directory to read benchmarks from: snakedata_r232 (Snakefile) or snakedata (Snakefile.r226) [default: %(default)s]')
parser.add_argument('--output-prefix', default=None,
                    help='write <prefix>.png and <prefix>.tsv [default: plots/singlem_benchmarks for snakedata, otherwise plots/singlem_benchmarks.<data dir>]')
args = parser.parse_args()
data_dir = args.data_dir.rstrip('/')
output_prefix = args.output_prefix or (
    'plots/singlem_benchmarks' if os.path.basename(data_dir) == 'snakedata'
    else 'plots/singlem_benchmarks.' + os.path.basename(data_dir))

rows = []
for stage in STAGES:
    for f in glob.glob(f'{data_dir}/benchmarks/{stage}/*/*.txt'):
        version = f.split('/')[-2]
        sample = os.path.basename(f)[:-4]
        b = pl.read_csv(f, separator='\t', null_values='NA')
        if b['max_rss'][0] is None:  # 0.19.0 has no long-read mode, so that job didn't really run
            continue
        rows.append(dict(stage=STAGES[stage], version=LABELS[version], sample=SAMPLES[sample],
                         wall_min=b['s'][0] / 60, rss_gb=b['max_rss'][0] / 1024))
if not rows:
    raise SystemExit(f'No benchmark files found under {data_dir}/benchmarks')
df = pl.DataFrame(rows)
df.sort(['stage', 'sample', 'version']).write_csv(output_prefix + '.tsv', separator='\t')

long = (df.unpivot(index=['stage', 'version', 'sample'], on=['wall_min', 'rss_gb'], variable_name='metric')
        .with_columns(pl.col('metric').replace({'wall_min': 'Wall time (minutes)', 'rss_gb': 'Peak RAM (GB)'}),
                      pl.col('value').map_elements(lambda v: f'{v:.2g}' if v < 10 else f'{v:.0f}', return_dtype=pl.String).alias('label')))
import pandas as pd
pdf = pd.DataFrame(long.to_dicts())
import pandas as pd
pdf['version'] = pd.Categorical(pdf['version'], [LABELS[v] for v in VERSIONS])
pdf['stage'] = pd.Categorical(pdf['stage'], list(STAGES.values()))
pdf['metric'] = pd.Categorical(pdf['metric'], ['Wall time (minutes)', 'Peak RAM (GB)'])
pdf['sample'] = pd.Categorical(pdf['sample'], list(SAMPLES.values()))

dodge = position_dodge(width=0.72)
pdf['dev_label'] = ['    ' + l if v == LABELS['dev'] else '' for l, v in zip(pdf['label'], pdf['version'])]
p = (ggplot(pdf, aes('sample', 'value', color='version', group='version'))
     + geom_point(position=dodge, size=3.2)
     + geom_text(aes(label='dev_label'), position=dodge, ha='left', va='center', size=7.5, color='#52514e')
     + facet_grid('metric ~ stage', scales='free', space='free_x')
     + scale_y_log10()
     + scale_x_discrete(expand=(0, 0.75))
     + scale_color_manual(values=COLORS)
     + labs(x='', y='', color='SingleM version',
            title='SingleM performance across versions (4 threads; log scale; labels = dev now)',
            caption=f'0.19.0 has no long-read mode, so it was not run on nanopore. Values for every version are in {output_prefix}.tsv ({data_dir})')
     + theme_minimal(base_size=9)
     + theme(figure_size=(12, 6.5), panel_grid_major_x=element_blank(), panel_grid_minor=element_blank(),
             panel_grid_major_y=element_line(color='#e6e5e0', size=0.4), legend_position='top',
             plot_background=element_rect(fill='#fcfcfb', color='#fcfcfb'),
             strip_text=element_text(weight='bold')))
p.save(output_prefix + '.png', dpi=150, verbose=False)
print(df.sort(['stage', 'sample', 'version']))
