#!/usr/bin/env python3
"""Regenerate AE traceability figures and tables added on 2026-06-09.
This script is intentionally lightweight: it reads the CSV files in ../data
and regenerates the traceability figures in ../figures.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
FIGS = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

heat = pd.read_csv(DATA/'ae_requirement_dataset_support_matrix_20260609.csv', index_col=0)
fig, ax = plt.subplots(figsize=(7.6, 4.4), dpi=300)
im = ax.imshow(heat.values, vmin=0, vmax=3, aspect='auto')
ax.set_xticks(range(len(heat.columns)), [c.replace(' ', '\n') if len(c)<12 else c for c in heat.columns], rotation=35, ha='right')
ax.set_yticks(range(len(heat.index)), heat.index)
ax.set_title('AE-requested data-to-parameter traceability map')
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        ax.text(j, i, f'{heat.values[i,j]:g}', ha='center', va='center', fontsize=7)
fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025, label='Support level: 0 none, 1 indirect, 2 component, 3 closest/direct')
fig.tight_layout()
fig.savefig(FIGS/'Figure S11 AE-requested traceability map.png', dpi=300)
plt.close(fig)

status = pd.read_csv(DATA/'traceability_validation_status_20260609.csv').sort_values('support_score_0_to_3')
fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=300)
y = np.arange(len(status))
ax.barh(y, status['support_score_0_to_3'])
ax.set_yticks(y, status['component'])
ax.set_xlim(0,3)
ax.set_xlabel('Evidence support score (0 = none, 3 = direct/closest)')
ax.set_title('Validation status by model component')
for yi, score, label in zip(y, status['support_score_0_to_3'], status['support_status']):
    ax.text(min(score+0.05,2.95), yi, label, va='center', fontsize=7)
fig.tight_layout()
fig.savefig(FIGS/'Figure S12 component validation status ladder.png', dpi=300)
plt.close(fig)

plt_data = pd.read_csv(DATA/'kulkarni2021_digitized_plate_load_data.csv')
fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=300)
for treatment, group in plt_data.groupby('treatment'):
    group = group.sort_values('settlement_mm')
    ax.plot(group['settlement_mm']/120.0, group['pressure_kpa'], marker='o', ms=3, label=treatment)
ax.set_xlabel('Normalized settlement, s/B (B = 120 mm)')
ax.set_ylabel('Pressure (kPa)')
ax.set_title('Digitized Kulkarni et al. (2021) PLT anchor')
ax.legend(title='State')
ax.grid(True, linewidth=0.3, alpha=0.5)
fig.tight_layout()
fig.savefig(FIGS/'Figure S13 normalized PLT load-settlement anchor.png', dpi=300)
plt.close(fig)

anchors = pd.read_csv(DATA/'degradation_retention_anchor_points_20260609.csv')
plot_df = anchors[(anchors['retention'] < 0.995) | (anchors['source'].str.contains('Ji|Tao', regex=True))].copy()
plot_df['label'] = plot_df['source'].str.replace(' et al.','', regex=False) + '\n' + plot_df['exposure_label'].str.replace(', ', '\n', regex=False)
plot_df = plot_df.sort_values('retention')
fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=300)
y = np.arange(len(plot_df))
ax.barh(y, plot_df['retention'])
ax.set_yticks(y, plot_df['label'])
ax.set_xlim(0,1.05)
ax.set_xlabel('Retained response ratio (dimensionless)')
ax.set_title('External degradation/retention anchors used as bounds')
for yi, val, metric in zip(y, plot_df['retention'], plot_df['metric']):
    ax.text(min(val+0.02,1.02), yi, f'{val:.2f} ({metric})', va='center', fontsize=7)
fig.tight_layout()
fig.savefig(FIGS/'Figure S14 degradation retention anchors.png', dpi=300)
plt.close(fig)
print('Traceability figures regenerated.')
